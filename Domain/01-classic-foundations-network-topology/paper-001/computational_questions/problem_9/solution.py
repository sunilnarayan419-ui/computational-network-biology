import numpy as np
import pandas as pd
import networkx as nx
from scipy import stats, optimize
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Configuration
INPUT = "edges.csv"              # columns: source,target[,time]
DIRECTED = False                 # True -> in-degree and out-degree reported separately
ANNOTATIONS = None               # optional CSV with a "node" column and annotation columns
EXCLUDE_NODES = []               # optional nodes to drop (e.g. currency metabolites)
TOP_N = 20
HUB_QUANTILE = 0.95              # candidate hubs: degree >= this quantile (ties included)
INITIAL_ATTRACTIVENESS = 0.0     # a in attachment rate (k + a)^alpha
BETWEENNESS_SAMPLES = 500        # sampled sources for betweenness on large networks
MAX_DRAW_NODES = 300
SEED = 42

# 1. Load and clean the edge list
edges = pd.read_csv(INPUT)
has_time = "time" in edges.columns
edges = edges[["source", "target"] + (["time"] if has_time else [])].dropna()
edges[["source", "target"]] = edges[["source", "target"]].astype(str)
edges = edges[(edges["source"] != edges["target"])
              & ~edges["source"].isin(EXCLUDE_NODES) & ~edges["target"].isin(EXCLUDE_NODES)]

if not DIRECTED:
    pairs = np.sort(edges[["source", "target"]].to_numpy(), axis=1)
    edges["source"], edges["target"] = pairs[:, 0], pairs[:, 1]

if has_time:
    edges["time"] = pd.to_numeric(edges["time"], errors="coerce")
    edges = edges.dropna(subset=["time"]).sort_values("time", kind="stable")

edges = edges.drop_duplicates(["source", "target"])     # earliest occurrence is kept

if edges.empty:
    raise ValueError("No valid edges after cleaning.")

G = nx.from_pandas_edgelist(edges, "source", "target",
                            create_using=nx.DiGraph if DIRECTED else nx.Graph)
print("Nodes:", G.number_of_nodes(), "| Edges:", G.number_of_edges(), "| Directed:", DIRECTED)

# 2. Degree, maximum-degree set (with ties), ranking and candidate hubs
if DIRECTED:
    deg = pd.DataFrame({"in_degree": dict(G.in_degree()), "out_degree": dict(G.out_degree())})
else:
    deg = pd.DataFrame({"degree": dict(G.degree())})

# 3. Other centralities, reported separately from degree
sampled = BETWEENNESS_SAMPLES if G.number_of_nodes() > BETWEENNESS_SAMPLES else None
deg["betweenness"] = pd.Series(nx.betweenness_centrality(G, k=sampled, seed=SEED))

try:
    deg["eigenvector"] = pd.Series(nx.eigenvector_centrality_numpy(G))
except Exception as err:
    deg["eigenvector"] = np.nan
    print("Eigenvector centrality not computed:", err)

if ANNOTATIONS:
    ann = pd.read_csv(ANNOTATIONS)
    ann["node"] = ann["node"].astype(str)
    deg = deg.join(ann.drop_duplicates("node").set_index("node"))

measures = ["in_degree", "out_degree"] if DIRECTED else ["degree"]

for m in measures:
    kmax = deg[m].max()
    deg[f"{m}_rank"] = deg[m].rank(method="min", ascending=False).astype(int)
    threshold = deg[m].quantile(HUB_QUANTILE)
    deg[f"{m}_hub"] = deg[m] >= threshold

    print(f"\n=== {m} ===")
    print(f"k_max = {kmax} | nodes sharing k_max: {', '.join(deg.index[deg[m] == kmax])}")
    print(f"Candidate hubs (>= {HUB_QUANTILE:.0%} quantile, k >= {threshold:g}): {int(deg[f'{m}_hub'].sum())}")
    print(deg.sort_values(m, ascending=False).head(TOP_N).round(4).to_string())

    rho, p = stats.spearmanr(deg[m], deg["betweenness"])
    print(f"Spearman {m} vs betweenness: rho = {rho:.3f} (degree is not interchangeable with centrality)")

    dist = deg[m].value_counts().sort_index().rename_axis(m).reset_index(name="nodes")
    dist["P(k)"] = dist["nodes"] / len(deg)
    dist.to_csv(f"degree_distribution_{m}.csv", index=False)

deg.sort_values(measures[0], ascending=False).to_csv("hub_table.csv")

# 4. Figure: degree distributions and, for small networks, the hub-highlighted network
draw = G.number_of_nodes() <= MAX_DRAW_NODES
fig, axes = plt.subplots(1, len(measures) + draw, figsize=(5 * (len(measures) + draw), 4), squeeze=False)

for ax, m in zip(axes[0], measures):
    counts = deg[m].value_counts().sort_index()
    ax.loglog(counts.index[counts.index > 0], counts[counts.index > 0] / len(deg), "o")
    ax.axvline(deg[m].quantile(HUB_QUANTILE), color="r", linestyle="--", label="hub threshold")
    ax.set(xlabel=m, ylabel="P(k)", title=f"{m} distribution")
    ax.legend()

if draw:
    ax = axes[0, -1]
    pos = nx.spring_layout(G, seed=SEED)
    hub = deg[f"{measures[0]}_hub"]
    nx.draw_networkx(G, pos, ax=ax, with_labels=False, node_size=20 + 20 * deg[measures[0]].reindex(list(G)).to_numpy(),
                     node_color=["tab:red" if hub[n] else "lightgray" for n in G], edge_color="lightgray")
    ax.set_title(f"Hubs highlighted ({measures[0]})")
    ax.axis("off")

plt.tight_layout()
plt.savefig("hubs.png", dpi=200)

# 5. Preferential attachment (requires time-stamped edges)
if not has_time:
    print("\nNo time column: preferential attachment cannot be tested from a static network.")
    print("Conclusions are limited to hub identification and topology description.")
else:
    # Attachment event = edge joining a node new to the network (source, if directed) to an existing node.
    # The chosen node's degree (in-degree if directed) is recorded before the edge's time step is applied.
    k_final = deg["in_degree" if DIRECTED else "degree"].max()
    count = np.zeros(int(k_final) + 2)                # count[k] = existing nodes with degree k
    deg_now, seen = {}, set()
    snapshots, n_events, sum_log, chosen_all = [], [], 0.0, []
    used = skipped = 0

    for _, g in edges.groupby("time", sort=True):
        chosen = []
        for u, v in zip(g["source"], g["target"]):
            if DIRECTED:
                old, new_ok = v, (u not in seen) and (v in seen)
            else:
                new_ok = (u in seen) != (v in seen)
                old = v if u not in seen else u
            if new_ok and (deg_now[old] > 0 or INITIAL_ATTRACTIVENESS > 0):
                chosen.append(deg_now[old])
            else:
                skipped += 1

        if chosen:
            snapshots.append(count.copy())
            n_events.append(len(chosen))
            sum_log += np.log(np.array(chosen) + INITIAL_ATTRACTIVENESS).sum()
            chosen_all += chosen
            used += len(chosen)

        for u, v in zip(g["source"], g["target"]):
            for node in (u, v):
                if node not in seen:
                    seen.add(node)
                    deg_now[node] = 0
                    count[0] += 1
            for node in ([v] if DIRECTED else [u, v]):
                count[deg_now[node]] -= 1
                deg_now[node] += 1
                count[deg_now[node]] += 1

    print(f"\nAttachment events used: {used} | edges not counted as events: {skipped}")

    if used < 30:
        print("Too few attachment events to estimate the attachment exponent.")
    else:
        Hm, mg = np.array(snapshots), np.array(n_events, float)
        ks = np.arange(Hm.shape[1]) + INITIAL_ATTRACTIVENESS

        def weights(alpha):
            with np.errstate(divide="ignore"):
                w = ks ** alpha
            if INITIAL_ATTRACTIVENESS == 0:
                w[0] = 0.0                    # nodes of degree 0 cannot be chosen when a = 0
            return w

        def loglik(alpha):
            return alpha * sum_log - (mg * np.log(Hm @ weights(alpha))).sum()

        # Maximum-likelihood attachment exponent: P(choose node i) = (k_i + a)^alpha / sum_j (k_j + a)^alpha
        lo, hi = -1.0, 4.0
        alpha_hat = optimize.minimize_scalar(lambda x: -loglik(x), bounds=(lo, hi), method="bounded").x
        ll_max = loglik(alpha_hat)

        def profile_bound(a, b):
            f = lambda x: loglik(x) - (ll_max - 1.92)       # 95% profile-likelihood interval
            return optimize.brentq(f, a, b) if f(a) * f(b) < 0 else np.nan

        ci = (profile_bound(lo, alpha_hat), profile_bound(alpha_hat, hi))
        print(f"alpha = {alpha_hat:.3f} (95% CI {ci[0]:.3f} to {ci[1]:.3f}), a = {INITIAL_ATTRACTIVENESS}")

        for name, a0 in [("uniform attachment (alpha = 0)", 0.0), ("linear preferential attachment (alpha = 1)", 1.0)]:
            lr = 2 * (ll_max - loglik(a0))
            print(f"LR test vs {name}: statistic = {lr:.2f}, p = {stats.chi2.sf(lr, 1):.4g}")

        # Observed vs model-predicted attachment counts by degree class
        bins = np.unique(np.r_[0, 2 ** np.arange(0, int(np.ceil(np.log2(Hm.shape[1]))) + 1)])
        klass = np.digitize(np.arange(Hm.shape[1]), bins)
        observed = np.bincount(chosen_all, minlength=Hm.shape[1])
        table = pd.DataFrame({"k_class": klass, "observed": observed})

        for label, a_ in [("uniform", 0.0), ("linear", 1.0), ("fitted", alpha_hat)]:
            w = weights(a_)
            table[label] = (mg[:, None] * Hm * w / (Hm @ w)[:, None]).sum(axis=0)

        summary = table.groupby("k_class").sum()
        summary.index = [f"k >= {bins[i - 1]}" for i in summary.index]
        summary.round(2).to_csv("attachment_pattern.csv")
        print(summary.round(1).to_string())
        print("Caveat: edge order is assumed to reflect growth order; missing or delayed edges bias alpha.")
