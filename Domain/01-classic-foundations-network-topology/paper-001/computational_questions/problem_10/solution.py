import time
from collections import Counter

import numpy as np
import pandas as pd
import networkx as nx
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Configuration
INPUT = "edges.csv"                  # columns: source,target
DIRECTED = True
MOTIFS = {                           # induced, weakly connected target patterns (all of size m, 3 <= m <= 6)
    "feed_forward_loop": [("a", "b"), ("b", "c"), ("a", "c")],
}
PILOT_SAMPLES = 2000                 # uniform samples used to estimate conditional probabilities
GUIDED_SAMPLES = 5000                # probability-guided samples used for frequency estimation
SMOOTHING = 2.0                      # Beta-prior strength (pseudo-observations) for unseen contexts
WEIGHT_FLOOR = 0.01                  # added to every extension weight so all configurations stay reachable
N_RANDOM = 20                        # randomized networks for the null distribution
SWAPS_PER_EDGE = 5
BOOTSTRAPS = 1000                    # bootstrap replicates for sampling uncertainty
SEED = 42

rng = np.random.default_rng(SEED)

# 1. Load network (self-loops removed, duplicate edges merged)
edges = pd.read_csv(INPUT, usecols=["source", "target"]).dropna().astype(str)
edges = edges[edges["source"] != edges["target"]]

if edges.empty:
    raise ValueError("No valid edges after removing missing values and self-loops.")

G = nx.from_pandas_edgelist(edges, "source", "target", create_using=nx.DiGraph if DIRECTED else nx.Graph)
print("Nodes:", G.number_of_nodes(), "| Edges:", G.number_of_edges(), "| Directed:", DIRECTED)
print("Weakly connected components:", nx.number_connected_components(G.to_undirected()))

# Targets are validated and stored as directed graphs (an undirected edge becomes two arcs)
targets = {}
for name, motif_edges in MOTIFS.items():
    T = nx.DiGraph(motif_edges) if DIRECTED else nx.Graph(motif_edges).to_directed()
    if not nx.is_weakly_connected(T) or not 3 <= T.number_of_nodes() <= 6:
        raise ValueError(f"Motif {name} must be connected with 3-6 nodes.")
    targets[name] = T

if len({T.number_of_nodes() for T in targets.values()}) > 1:
    raise ValueError("All target motifs must have the same size.")


def expand(D, S):
    """Frontier of node set S: {u: [arcs u->S, arcs S->u]} and number of arcs inside S."""
    frontier, inside = {}, 0
    for s in S:
        for u in D.predecessors(s):
            if u in S:
                inside += 1
            else:
                frontier.setdefault(u, [0, 0])[0] += 1
        for u in D.successors(s):
            if u not in S:
                frontier.setdefault(u, [0, 0])[1] += 1
    return frontier, inside


def draw(D, nodes, m, weight):
    """Grow a connected m-node set from a uniform seed; next node chosen with probability proportional to weight."""
    S = {nodes[rng.integers(len(nodes))]}
    path = []

    for _ in range(m - 1):
        frontier, inside = expand(D, S)
        if not frontier:
            return None
        candidates = list(frontier)
        contexts = [(len(S), inside, *frontier[u]) for u in candidates]
        w = np.array([weight(c) for c in contexts])
        i = rng.choice(len(candidates), p=w / w.sum())
        path.append(contexts[i])
        S.add(candidates[i])

    return S, path


def set_probability(D, X, n_seed, weight):
    """Exact probability that draw() produces node set X, summed over all valid construction orders."""
    X, m = list(X), len(X)
    level = {1 << i: 1.0 / n_seed for i in range(m)}

    for size in range(1, m):
        nxt = {}
        for mask, p_mask in level.items():
            S = {X[i] for i in range(m) if mask >> i & 1}
            frontier, inside = expand(D, S)
            total = sum(weight((size, inside, *ab)) for ab in frontier.values())

            for j in range(m):
                if mask >> j & 1 or X[j] not in frontier:
                    continue
                step = weight((size, inside, *frontier[X[j]])) / total
                nxt[mask | 1 << j] = nxt.get(mask | 1 << j, 0.0) + p_mask * step

        level = nxt

    return level.get((1 << m) - 1, 0.0)


def search(D, T):
    """Pilot -> conditional probabilities -> guided importance sampling of the target motif."""
    nodes = [v for v in D if D.degree(v) > 0]
    m = T.number_of_nodes()

    # Stage A: uniform pilot. Context B = (|S|, arcs inside S, arcs u->S, arcs S->u); A = final set is the motif.
    n_ctx, hit_ctx, n_ok, n_hit = Counter(), Counter(), 0, 0
    for _ in range(PILOT_SAMPLES):
        out = draw(D, nodes, m, lambda c: 1.0)
        if out is None:
            continue
        X, path = out
        hit = nx.is_isomorphic(D.subgraph(X), T)
        n_ok += 1
        n_hit += hit
        for c in path:
            n_ctx[c] += 1
            hit_ctx[c] += hit

    p0 = n_hit / n_ok if n_ok else 0.0
    guided = p0 > 0
    if guided:
        # Smoothed conditional probability P(motif | context); unseen contexts fall back to p0
        weight = lambda c: (hit_ctx[c] + SMOOTHING * p0) / (n_ctx[c] + SMOOTHING) + WEIGHT_FLOOR
    else:
        weight = lambda c: 1.0

    # Stage B: guided sampling; 1/p(X) corrects for unequal sampling probabilities
    inv_p, is_hit, found = np.zeros(GUIDED_SAMPLES), np.zeros(GUIDED_SAMPLES), set()
    for i in range(GUIDED_SAMPLES):
        out = draw(D, nodes, m, weight)
        if out is None:
            continue
        X = out[0]
        inv_p[i] = 1.0 / set_probability(D, X, len(nodes), weight)
        if nx.is_isomorphic(D.subgraph(X), T):
            is_hit[i] = 1
            found.add(frozenset(X))

    num, den = inv_p * is_hit, inv_p
    count_hat, total_hat = num.mean(), den.mean()
    freq_hat = count_hat / total_hat if total_hat > 0 else np.nan

    # Bootstrap over sampled configurations (independent draws) for sampling uncertainty
    idx = rng.integers(GUIDED_SAMPLES, size=(BOOTSTRAPS, GUIDED_SAMPLES))
    b_num, b_den = num[idx].mean(axis=1), den[idx].mean(axis=1)
    with np.errstate(all="ignore"):
        b_freq = b_num / b_den

    table = pd.DataFrame([
        {"context (|S|, arcs in S, u->S, S->u)": c, "n_B": n_ctx[c], "n_AB": hit_ctx[c],
         "P_raw": hit_ctx[c] / n_ctx[c],
         "P_smoothed": (hit_ctx[c] + SMOOTHING * p0) / (n_ctx[c] + SMOOTHING),
         "ci_low": stats.beta.ppf(0.025, hit_ctx[c] + 0.5, n_ctx[c] - hit_ctx[c] + 0.5),
         "ci_high": stats.beta.ppf(0.975, hit_ctx[c] + 0.5, n_ctx[c] - hit_ctx[c] + 0.5)}
        for c in n_ctx
    ]).sort_values("P_smoothed", ascending=False) if n_ctx else pd.DataFrame()

    weights_ok = inv_p[inv_p > 0]
    return {
        "count": count_hat, "count_ci": np.percentile(b_num, [2.5, 97.5]),
        "total": total_hat, "freq": freq_hat, "freq_ci": np.nanpercentile(b_freq, [2.5, 97.5]),
        "ess": weights_ok.sum() ** 2 / (weights_ok ** 2).sum() if weights_ok.size else 0.0,
        "pilot_hits": n_hit, "pilot_ok": n_ok, "guided": guided, "table": table, "found": found,
    }


# 2. Observed network and randomized (degree-preserving) null networks
D_obs = G if DIRECTED else G.to_directed()
results, null_counts = {}, {name: [] for name in targets}

for name, T in targets.items():
    start = time.perf_counter()
    results[name] = search(D_obs, T)
    results[name]["seconds"] = time.perf_counter() - start

n_swap = SWAPS_PER_EDGE * G.number_of_edges()
for _ in range(N_RANDOM):
    Gr = G.copy()
    seed = int(rng.integers(1_000_000_000))
    if DIRECTED:
        nx.directed_edge_swap(Gr, nswap=n_swap, max_tries=100 * n_swap, seed=seed)
    else:
        nx.double_edge_swap(Gr, nswap=n_swap, max_tries=100 * n_swap, seed=seed)

    D_rnd = Gr if DIRECTED else Gr.to_directed()
    for name, T in targets.items():
        null_counts[name].append(search(D_rnd, T)["count"])        # same pipeline as for the observed network

# 3. Report, significance, multiple-testing correction
p_values = []
for name in targets:
    obs, null = results[name]["count"], np.array(null_counts[name])
    p_values.append((1 + (null >= obs).sum()) / (1 + null.size))

adjusted = multipletests(p_values, method="fdr_bh")[1] if len(p_values) > 1 else p_values

print(f"\nPilot samples: {PILOT_SAMPLES} | guided samples: {GUIDED_SAMPLES} | null networks: {N_RANDOM} | seed {SEED}")

for (name, res), p, q in zip(results.items(), p_values, adjusted):
    null = np.array(null_counts[name])
    print(f"\n=== {name} (m = {targets[name].number_of_nodes()}) ===")
    print(f"Pilot: {res['pilot_hits']} motif hits in {res['pilot_ok']} uniform samples"
          + ("" if res["guided"] else " -> no hits, guided stage fell back to uniform sampling"))
    print(res["table"].head(10).round(4).to_string(index=False))
    res["table"].to_csv(f"conditional_probabilities_{name}.csv", index=False)

    pd.DataFrame({"candidate_node_set": [";".join(sorted(map(str, s))) for s in res["found"]]}) \
        .to_csv(f"candidate_motifs_{name}.csv", index=False)

    print(f"Distinct matching node sets sampled (candidate discovery): {len(res['found'])}")
    print(f"Estimated motif count = {res['count']:.1f} (95% bootstrap CI {res['count_ci'][0]:.1f} to {res['count_ci'][1]:.1f})")
    print(f"Estimated frequency among connected {targets[name].number_of_nodes()}-node sets = {res['freq']:.3g} "
          f"(95% CI {res['freq_ci'][0]:.3g} to {res['freq_ci'][1]:.3g}); effective sample size = {res['ess']:.1f}")
    print(f"Null estimate: mean = {null.mean():.1f}, sd = {null.std(ddof=1) if null.size > 1 else np.nan:.1f}; "
          f"empirical p = {p:.4f}, BH-adjusted p = {q:.4f}")
    print(f"Runtime of search on observed network: {res['seconds']:.1f} s")

print("\nNote: the smallest attainable empirical p-value is 1/(N_RANDOM + 1); use more null networks for strong claims.")
print("Conditional-probability intervals treat pilot samples as independent although neighbourhoods overlap.")

# 4. Observed estimate vs null distribution
fig, axes = plt.subplots(1, len(targets), figsize=(5 * len(targets), 4), squeeze=False)
for ax, name in zip(axes[0], targets):
    ax.hist(null_counts[name], bins=10, color="lightgray", label="randomized networks")
    ax.axvline(results[name]["count"], color="r", label="observed estimate")
    ax.axvspan(*results[name]["count_ci"], color="r", alpha=0.15)
    ax.set(xlabel="estimated motif count", ylabel="networks", title=name)
    ax.legend(fontsize=7)
plt.tight_layout()
plt.savefig("motif_null_comparison.png", dpi=200)
