from collections import Counter
from itertools import combinations

import cobra
import numpy as np
import pandas as pd
import networkx as nx
from sklearn.metrics import normalized_mutual_info_score as nmi
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Configuration
MODEL = "model.xml"            # SBML genome-scale reconstruction of a cellular metabolic network
CURRENCY = []                  # metabolite IDs removed before projection (sensitivity run vs. all metabolites)
EXCLUDE_REACTIONS = []         # e.g. biomass reaction IDs; boundary (exchange/demand/sink) reactions are always excluded
RESOLUTIONS = [0.5, 1.0, 2.0, 4.0]   # Louvain resolution, coarse -> fine
N_RUNS = 3                     # repeated Louvain runs per resolution (stability)
N_RANDOM = 10                  # randomized metabolic networks
SWAPS_PER_EDGE = 5
NESTED_MIN = 0.80              # screening thresholds for the signature summary
STABILITY_MIN = 0.70
SEED = 42

rng = np.random.default_rng(SEED)

# 1. Reaction-metabolite incidence from the stoichiometric model (negative = substrate, positive = product)
model = cobra.io.read_sbml_model(MODEL)
subsystem = {r.id: r.subsystem for r in model.reactions if r.subsystem}
rows = [(r.id, m.id, c) for r in model.reactions
        if not r.boundary and r.id not in EXCLUDE_REACTIONS for m, c in r.metabolites.items()]
inc_all = pd.DataFrame(rows, columns=["reaction", "metabolite", "coefficient"])

if inc_all.empty:
    raise ValueError("No internal reactions found in the model.")

inc = inc_all[~inc_all["metabolite"].isin(CURRENCY)]
inc.to_csv("bipartite_incidence.csv", index=False)

print(f"Model {model.id}: {inc_all['reaction'].nunique()} internal reactions, "
      f"{inc_all['metabolite'].nunique()} metabolites, compartments: {sorted(model.compartments)}")
print(f"Reactions with a subsystem annotation: {len(set(inc_all['reaction']) & set(subsystem))}")
print(f"Currency metabolites removed: {len(CURRENCY)} | incidences analysed: {len(inc)}")


def project_reactions(table):
    """Reaction graph; edge weight = number of metabolites shared by two reactions (any role)."""
    W = Counter()
    for _, g in table.groupby("metabolite"):
        W.update(combinations(sorted(g["reaction"].unique()), 2))
    G = nx.Graph()
    G.add_weighted_edges_from((a, b, w) for (a, b), w in W.items())
    return G


def project_metabolites(table):
    """Metabolite graph; edge weight = number of reactions converting one metabolite into the other."""
    W = Counter()
    for _, g in table.groupby("reaction"):
        for s in g.loc[g["coefficient"] < 0, "metabolite"]:
            for p in g.loc[g["coefficient"] > 0, "metabolite"]:
                if s != p:
                    W[tuple(sorted((s, p)))] += 1
    G = nx.Graph()
    G.add_weighted_edges_from((a, b, w) for (a, b), w in W.items())
    return G


def bipartite_null(table):
    """Degree-preserving swaps of reaction-metabolite incidences; each reaction keeps its substrate/product roles."""
    r = table["reaction"].to_numpy()
    m = table["metabolite"].to_numpy().copy()
    pairs = set(zip(r, m))

    for _ in range(SWAPS_PER_EDGE * len(r)):
        i, j = rng.integers(len(r), size=2)
        if r[i] == r[j] or m[i] == m[j] or (r[i], m[j]) in pairs or (r[j], m[i]) in pairs:
            continue
        pairs -= {(r[i], m[i]), (r[j], m[j])}
        pairs |= {(r[i], m[j]), (r[j], m[i])}
        m[i], m[j] = m[j], m[i]

    return pd.DataFrame({"reaction": r, "metabolite": m, "coefficient": table["coefficient"].to_numpy()})


def pipeline(table):
    """Multi-resolution Louvain modules and their modularity, stability, nestedness and annotation agreement."""
    stats_out, partitions, graphs = {}, {}, {}

    for rep, G in (("reaction", project_reactions(table)), ("metabolite", project_metabolites(table))):
        if G.number_of_edges() == 0:
            continue
        graphs[rep], nodes, previous = G, list(G), None

        for res in RESOLUTIONS:
            runs = [nx.community.louvain_communities(G, weight="weight", resolution=res, seed=SEED + i)
                    for i in range(N_RUNS)]
            labels = [[next(k for k, c in enumerate(comms) if n in c) for n in nodes] for comms in runs]
            current = dict(zip(nodes, labels[0]))
            partitions[(rep, res)] = current

            stats_out[(rep, "Q", res)] = nx.community.modularity(G, runs[0], weight="weight")
            if N_RUNS > 1:
                stats_out[(rep, "stability", res)] = np.mean([nmi(labels[a], labels[b])
                                                              for a, b in combinations(range(N_RUNS), 2)])
            if previous is not None:            # share of fine-module nodes inside their best-matching coarse module
                ct = pd.crosstab(pd.Series(current), pd.Series(previous))
                stats_out[(rep, "nested", res)] = ct.max(axis=1).sum() / ct.to_numpy().sum()
            if rep == "reaction":               # agreement with curated subsystems (annotated reactions only)
                annotated = [n for n in nodes if n in subsystem]
                if len(annotated) > 1:
                    stats_out[(rep, "nmi_subsystem", res)] = nmi([subsystem[n] for n in annotated],
                                                                 [current[n] for n in annotated])
            previous = current

    return stats_out, partitions, graphs


# 2. Observed analysis
observed, partitions, graphs = pipeline(inc)

if not graphs:
    raise ValueError("Projection produced no edges; check the model and CURRENCY list.")

conn_rows, parent_rows, assign = [], [], {}
for (rep, res), part in partitions.items():
    G = graphs[rep]
    e = pd.DataFrame([(part[u], part[v], d["weight"]) for u, v, d in G.edges(data=True)], columns=["cu", "cv", "w"])
    internal = e[e["cu"] == e["cv"]].groupby("cu")["w"].sum()
    cross = e[e["cu"] != e["cv"]]
    external = pd.concat([cross.groupby("cu")["w"].sum(), cross.groupby("cv")["w"].sum()]).groupby(level=0).sum()
    table = pd.DataFrame({"size": pd.Series(part).value_counts(), "internal": internal,
                          "external": external}).fillna(0)
    table["internal_fraction"] = table["internal"] / (table["internal"] + table["external"])
    table.insert(0, "resolution", res)
    table.insert(0, "representation", rep)
    conn_rows.append(table.reset_index(names="module"))

    assign.setdefault(rep, pd.DataFrame(index=list(G)))[f"module_res_{res}"] = pd.Series(part)

    prev_res = RESOLUTIONS[RESOLUTIONS.index(res) - 1] if res != RESOLUTIONS[0] else None
    if prev_res is not None:
        ct = pd.crosstab(pd.Series(part), pd.Series(partitions[(rep, prev_res)]))
        parent_rows.append(pd.DataFrame({"representation": rep, "coarse_resolution": prev_res, "fine_resolution": res,
                                         "fine_module": ct.index, "parent_module": ct.idxmax(axis=1).to_numpy(),
                                         "fraction_in_parent": (ct.max(axis=1) / ct.sum(axis=1)).to_numpy()}))

connectivity = pd.concat(conn_rows, ignore_index=True)
connectivity.to_csv("module_connectivity.csv", index=False)
pd.concat(parent_rows, ignore_index=True).to_csv("module_parent_child.csv", index=False)

for rep, frame in assign.items():
    if rep == "reaction":
        frame["subsystem"] = [subsystem.get(n) for n in frame.index]
    frame.to_csv(f"module_assignments_{rep}.csv")

# 3. Randomized metabolic networks (same pipeline)
null = {k: [] for k in observed}
for _ in range(N_RANDOM):
    null_stats, _, _ = pipeline(bipartite_null(inc))
    for k in observed:
        null[k].append(null_stats.get(k, np.nan))

# 4. Report
report = pd.DataFrame([
    {"representation": k[0], "metric": k[1], "resolution": k[2], "observed": v,
     "null_mean": np.nanmean(null[k]), "null_sd": np.nanstd(null[k], ddof=1) if N_RANDOM > 1 else np.nan,
     "empirical_p": (1 + np.nansum(np.array(null[k]) >= v)) / (1 + N_RANDOM)}
    for k, v in observed.items()
])
report["z"] = (report["observed"] - report["null_mean"]) / report["null_sd"]
report.to_csv("hierarchical_modularity_report.csv", index=False)

print(f"\nLouvain resolutions {RESOLUTIONS} | {N_RUNS} runs each | {N_RANDOM} randomized networks | seed {SEED}")
print(report.round(4).to_string(index=False))

summary = connectivity.groupby(["representation", "resolution"]).agg(
    modules=("module", "count"), median_size=("size", "median"), mean_internal_fraction=("internal_fraction", "mean"))
print("\nModule structure by resolution:\n", summary.round(3).to_string())

for rep in graphs:
    r = report[report["representation"] == rep]
    q, nest, stab, ann = (r[r["metric"] == m] for m in ("Q", "nested", "stability", "nmi_subsystem"))
    signatures = {
        "modularity above null at every resolution": bool((q["empirical_p"] <= 0.05).all()),
        f"nestedness >= {NESTED_MIN} and above null mean": bool(len(nest) and (nest["observed"] >= NESTED_MIN).all()
                                                               and (nest["observed"] > nest["null_mean"]).all()),
        f"stability >= {STABILITY_MIN} at every resolution": bool(len(stab) and (stab["observed"] >= STABILITY_MIN).all()),
    }
    if rep == "reaction" and len(ann):
        signatures["subsystem agreement above null at >= 1 resolution"] = bool((ann["empirical_p"] <= 0.05).any())

    print(f"\nSignatures for {rep} network ({sum(signatures.values())}/{len(signatures)} met):")
    for name, ok in signatures.items():
        print(f"  {'yes' if ok else 'no '}  {name}")

print("\nThese are screening criteria. Modularity at one resolution alone does not demonstrate hierarchical modularity.")

# 5. Optional currency-metabolite sensitivity (observed network only)
if CURRENCY:
    alt, _, _ = pipeline(inc_all)
    sens = pd.DataFrame([{"representation": k[0], "metric": k[1], "resolution": k[2],
                          "currency_removed": observed[k], "all_metabolites": alt.get(k, np.nan)} for k in observed])
    print("\nCurrency-metabolite sensitivity:\n", sens.round(4).to_string(index=False))

# 6. Figure: modularity and nestedness, observed vs randomized
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for rep in graphs:
    for ax, metric in zip(axes, ("Q", "nested")):
        r = report[(report["representation"] == rep) & (report["metric"] == metric)]
        if r.empty:
            continue
        ax.plot(r["resolution"], r["observed"], "o-", label=f"{rep} observed")
        ax.errorbar(r["resolution"], r["null_mean"], r["null_sd"], fmt="s--", capsize=3, label=f"{rep} randomized")
for ax, label in zip(axes, ("modularity Q", "nestedness between consecutive resolutions")):
    ax.set(xscale="log", xlabel="Louvain resolution", ylabel=label)
    ax.legend(fontsize=7)
plt.tight_layout()
plt.savefig("hierarchical_modularity.png", dpi=200)
