import re
from collections import Counter, defaultdict

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import powerlaw

# Configuration
INPUT = "reactions.csv"      # columns: reaction_id, substrates, products [, reversible]
SEP = ";"                    # separator for multiple metabolites inside one cell
PYRUVATE_IDS = ["pyruvate", "pyr_c"]   # must match identifiers used in INPUT (case-insensitive)
EXCLUDE_CURRENCY = True
CURRENCY = {"h2o", "h", "atp", "adp", "pi", "nad", "nadh", "nadp", "nadph"}
N_HUBS = 15
N_RANDOM = 100               # degree-preserving randomizations
SWAPS_PER_EDGE = 10
BETWEENNESS_K = None         # None = exact; integer = sampled pivots for large networks
SEED = 42
OUTPUT_TABLE = "metabolite_metrics.csv"
OUTPUT_FIG = "metabolic_topology.png"

# 1. Load reaction data
df = pd.read_csv(INPUT, dtype=str).fillna("")
if not {"reaction_id", "substrates", "products"} <= set(df.columns):
    raise ValueError("Required columns: reaction_id, substrates, products.")
if "reversible" not in df:
    df["reversible"] = "false"   # assumption when directionality is not supplied
n_input = len(df)
df = df.drop_duplicates("reaction_id")
print(f"Reactions: {len(df)} (duplicate reaction IDs removed: {n_input - len(df)})")


def split(cell):
    items = [x.strip() for x in cell.split(SEP) if x.strip()]
    return [x for x in items if not (EXCLUDE_CURRENCY and re.sub(r"_[a-z]$", "", x.lower()) in CURRENCY)]


# 2. Build the metabolite network: every substrate -> every product of the same reaction
D = nx.DiGraph()
rxns, produced, consumed = defaultdict(set), Counter(), Counter()
for r in df.itertuples():
    S, P = split(r.substrates), split(r.products)
    rev = r.reversible.strip().lower() in {"true", "1", "yes", "y"}
    D.add_nodes_from(S + P)
    for x in set(S + P):
        rxns[x].add(r.reaction_id)
    produced.update(set(P) | (set(S) if rev else set()))
    consumed.update(set(S) | (set(P) if rev else set()))
    for s in S:
        for p in P:
            if s != p:
                D.add_edge(s, p)
                if rev:
                    D.add_edge(p, s)
if D.number_of_nodes() == 0:
    raise ValueError("No metabolites found.")
U = D.to_undirected()

# 3. Connectivity and centrality
m = pd.DataFrame({
    "degree": dict(U.degree()),                 # distinct neighbouring metabolites
    "in_degree": dict(D.in_degree()),           # distinct metabolites converted into this one
    "out_degree": dict(D.out_degree()),         # distinct metabolites this one is converted into
    "n_reactions": {k: len(v) for k, v in rxns.items()},
    "produced_by": dict(produced),              # reactions producing it (reversible counted both ways)
    "consumed_by": dict(consumed),
}).fillna(0)
cnt = ["degree", "in_degree", "out_degree", "n_reactions", "produced_by", "consumed_by"]
m[cnt] = m[cnt].astype(int)
m["degree_centrality"] = pd.Series(nx.degree_centrality(U))
m["betweenness"] = pd.Series(nx.betweenness_centrality(U, k=BETWEENNESS_K, seed=SEED))
m = m.sort_values(["degree", "n_reactions"], ascending=False)
m.to_csv(OUTPUT_TABLE, index_label="metabolite")

print(f"\nMetabolites: {U.number_of_nodes()} | Neighbour links: {U.number_of_edges()}")
print(f"\nTop {N_HUBS} metabolites by degree:")
print(m.head(N_HUBS).round(4).to_string())

# 4. Pyruvate analysis (calculated, not assumed)
lookup = {x.lower(): x for x in m.index}
pyr_nodes = [lookup[p.lower()] for p in PYRUVATE_IDS if p.lower() in lookup]
rank = m["degree"].rank(ascending=False, method="min")
for x in pyr_nodes:
    print(f"\n{x}: degree {m.at[x, 'degree']} (rank {int(rank[x])} of {len(m)}, "
          f"above {(m['degree'] < m.at[x, 'degree']).mean() * 100:.1f}% of metabolites)")
    print(m.loc[x].round(4).to_dict())
    print("Reactions:", sorted(rxns[x]))
    print("Neighbours:", sorted(U[x]))
if not pyr_nodes:
    print("\nPyruvate identifier not found; set PYRUVATE_IDS.")

# 5. Degree distribution P(k) = N(k) / N and power-law analysis
def pk(values):
    c = pd.Series(values)[lambda s: s > 0].value_counts().sort_index()
    return c / len(values)


deg = m.loc[m["degree"] > 0, "degree"].to_numpy()
print("\nPower-law analysis (discrete MLE, xmin by KS minimisation):")
if len(set(deg)) >= 5:
    fit = powerlaw.Fit(deg, discrete=True, verbose=False)
    print(f"alpha = {fit.power_law.alpha:.3f} +/- {fit.power_law.sigma:.3f}, xmin = {fit.xmin}, "
          f"KS distance = {fit.power_law.D:.3f}, tail size = {int((deg >= fit.xmin).sum())} of {len(deg)}")
    for alt in ("exponential", "lognormal", "truncated_power_law"):
        try:
            R, p = fit.distribution_compare("power_law", alt, normalized_ratio=True)
            print(f"power_law vs {alt}: R = {R:.3f}, p = {p:.3f} (R > 0 favours power law)")
        except Exception as e:
            print(f"power_law vs {alt}: comparison failed ({e})")
    print("Note: no goodness-of-fit bootstrap is run; small tails give weak evidence either way.")
else:
    print("Too few distinct degree values for a power-law fit.")

# 6. Degree-preserving randomization
funcs = {"avg_clustering": nx.average_clustering, "degree_assortativity": nx.degree_assortativity_coefficient}
n_edges = U.number_of_edges()
if n_edges >= 10 and U.number_of_nodes() >= 4:
    null = {k: [] for k in funcs}
    for i in range(N_RANDOM):
        R = U.copy()
        nx.double_edge_swap(R, SWAPS_PER_EDGE * n_edges, max_tries=1000 * SWAPS_PER_EDGE * n_edges, seed=SEED + i)
        for k, f in funcs.items():
            null[k].append(f(R))
    print(f"\nObserved vs {N_RANDOM} degree-preserving randomizations:")
    for k, f in funcs.items():
        o, r = f(U), np.array(null[k])[np.isfinite(null[k])]
        p = (1 + np.sum(np.abs(r - r.mean()) >= abs(o - r.mean()))) / (len(r) + 1) if len(r) else np.nan
        print(f"{k}: observed {o:.4f}, random mean {r.mean():.4f} (sd {r.std():.4f}), empirical p = {p:.3f}")
else:
    print("\nNetwork too small for degree-preserving randomization.")

# 7. Visualization: P(k) and hub sub-network
fig, ax = plt.subplots(1, 2, figsize=(14, 6))
for col, style in (("degree", "o"), ("in_degree", "s"), ("out_degree", "^")):
    p = pk(m[col].to_numpy())
    ax[0].loglog(p.index, p.values, style, alpha=0.7, label=col)
ax[0].set(xlabel="k", ylabel="P(k)", title="Degree distributions")
ax[0].legend()
top = list(m.index[:N_HUBS]) + [x for x in pyr_nodes if x not in m.index[:N_HUBS]]
H = U.subgraph(top)
nx.draw_networkx(H, pos=nx.spring_layout(H, seed=SEED), ax=ax[1], font_size=7, width=0.4,
                 node_color=["gold" if n in pyr_nodes else "tab:red" for n in H],
                 node_size=[40 + 8 * U.degree(n) for n in H])
ax[1].set_title(f"Top {N_HUBS} hubs (pyruvate in gold)")
ax[1].axis("off")
fig.tight_layout()
fig.savefig(OUTPUT_FIG, dpi=150)
print(f"\nSaved {OUTPUT_TABLE} and {OUTPUT_FIG}")
