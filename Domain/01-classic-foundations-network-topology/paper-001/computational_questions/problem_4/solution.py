import re
from collections import defaultdict
from itertools import combinations

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

# Configuration
INPUT = "reactions.csv"      # columns: reaction_id, substrates, products [, reversible, pathway]
SEP = ";"
PYRUVATE_IDS = ["pyruvate", "pyr_c"]   # must match identifiers used in INPUT (case-insensitive)
CURRENCY = {"h2o", "h", "atp", "adp", "pi", "nad", "nadh", "nadp", "nadph"}  # never create links
PYRUVATE_LINKS = True        # False: pyruvate itself does not link reactions
DEPTH = 2                    # reaction-network steps around pyruvate-containing seed reactions
N_RANDOM = 200               # Erdos-Renyi realizations
SEED = 42
OUTPUT_FIG = "pyruvate_reaction_network_vs_er.png"


def base(x):
    return re.sub(r"_[a-z]$", "", x.lower())


# 1. Load reactions and collect the metabolites of each reaction
df = pd.read_csv(INPUT, dtype=str).fillna("")
if not {"reaction_id", "substrates", "products"} <= set(df.columns):
    raise ValueError("Required columns: reaction_id, substrates, products.")
n_input = len(df)
df = df.drop_duplicates("reaction_id")
print(f"Reactions: {len(df)} (duplicate reaction IDs removed: {n_input - len(df)})")

rxn_mets = {r.reaction_id: {x.strip() for cell in (r.substrates, r.products)
                            for x in cell.split(SEP) if x.strip()} for r in df.itertuples()}
pyr = {x for ms in rxn_mets.values() for x in ms if x.lower() in {p.lower() for p in PYRUVATE_IDS}}
seeds = {r for r, ms in rxn_mets.items() if ms & pyr}
if not seeds:
    raise ValueError("No pyruvate-containing reactions found; set PYRUVATE_IDS.")

# 2. Reaction-node projection: two reactions are linked if they share a non-currency metabolite
linkers = defaultdict(set)
for r, ms in rxn_mets.items():
    for x in ms:
        if base(x) not in CURRENCY and (PYRUVATE_LINKS or x not in pyr):
            linkers[x].add(r)
full = nx.Graph()
full.add_nodes_from(rxn_mets)
for x, rs in linkers.items():
    for a, b in combinations(sorted(rs), 2):
        if full.has_edge(a, b):
            full[a][b]["shared"].add(x)
        else:
            full.add_edge(a, b, shared={x})

# Pyruvate-associated reactions: seeds plus reactions within DEPTH steps of them
reach = nx.multi_source_dijkstra_path_length(full, seeds, cutoff=DEPTH)
G = full.subgraph(reach).copy()
N, L = G.number_of_nodes(), G.number_of_edges()
if N < 3:
    raise ValueError("Fewer than 3 reaction nodes; random-network comparison is not meaningful.")
p = 2 * L / (N * (N - 1))
print(f"Seed reactions containing pyruvate: {len(seeds)} | reaction nodes: {N} | edges: {L}")
print(f"Erdos-Renyi probability p = 2L/[N(N-1)] = {p:.4f}; expected edges = {p * N * (N - 1) / 2:.1f}")


# 3. Erdos-Renyi realizations and network statistics
def stats(g):
    comps = list(nx.connected_components(g))
    lcc = g.subgraph(max(comps, key=len))
    deg = np.array([d for _, d in g.degree()])
    return {"n_edges": g.number_of_edges(), "mean_degree": deg.mean(), "max_degree": deg.max(),
            "avg_clustering": nx.average_clustering(g), "n_components": len(comps),
            "largest_component_size": len(lcc), "isolated_nodes": sum(1 for _ in nx.isolates(g)),
            "avg_path_length_largest_component":
                nx.average_shortest_path_length(lcc) if len(lcc) > 1 else np.nan}


graphs = [nx.gnp_random_graph(N, p, seed=SEED + i) for i in range(N_RANDOM)]
obs = stats(G)
rand = pd.DataFrame([stats(g) for g in graphs])

rows = []
for k, o in obs.items():
    r = rand[k].dropna().to_numpy()
    p_emp = (1 + np.sum(np.abs(r - r.mean()) >= abs(o - r.mean()))) / (len(r) + 1) if len(r) and not np.isnan(o) else np.nan
    rows.append({"metric": k, "reaction_network": o, "ER_mean": r.mean() if len(r) else np.nan,
                 "ER_sd": r.std() if len(r) else np.nan, "empirical_p": p_emp})
print(f"\nComparison with {N_RANDOM} Erdos-Renyi realizations (empirical two-sided p):")
print(pd.DataFrame(rows).round(4).to_string(index=False))
print("Path lengths are computed within each graph's largest connected component.")


# 4. Degree distributions P(k) = N(k) / N
def pk(g):
    return pd.Series([d for _, d in g.degree()]).value_counts(normalize=True).sort_index()


pk_obs = pk(G)
pk_rand = pd.concat([pk(g) for g in graphs], axis=1).fillna(0).mean(axis=1)

# 5. Visualization
fig, ax = plt.subplots(1, 3, figsize=(16, 5))
for a, g, title in ((ax[0], G, "Pyruvate reaction-node projection"),
                    (ax[1], graphs[0], f"Erdos-Renyi G(N={N}, p={p:.3f})")):
    nx.draw_networkx(g, pos=nx.spring_layout(g, seed=SEED), ax=a, node_size=40, width=0.5, with_labels=False)
    a.set_title(title)
    a.axis("off")
ax[2].plot(pk_obs.index, pk_obs.values, "o-", label="Reaction network")
ax[2].plot(pk_rand.index, pk_rand.values, "s--", label=f"ER mean ({N_RANDOM})")
ax[2].set(xlabel="k", ylabel="P(k)", title="Degree distribution")
ax[2].legend()
fig.tight_layout()
fig.savefig(OUTPUT_FIG, dpi=150)
print("Saved", OUTPUT_FIG)
