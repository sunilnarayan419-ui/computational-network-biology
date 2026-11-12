from itertools import combinations

import numpy as np
import pandas as pd
import networkx as nx
from scipy import stats
from scipy.sparse.csgraph import shortest_path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Configuration
INPUT = "reactions.csv"            # columns: species,reaction,metabolite,coefficient[,group]
PROJECTION = "shared_metabolite"   # "shared_metabolite" or "substrate_product"
CURRENCY = []                      # metabolite IDs for the currency-handling sensitivity run
NULL_MODEL = "bipartite"           # "bipartite" (swap incidences, re-project) or "projection"
N_RANDOM = 20                      # randomized reference networks per species
SWAPS_PER_EDGE = 5
SEED = 42

rng = np.random.default_rng(SEED)

# 1. Load reaction-metabolite incidence data (negative coefficient = substrate, positive = product)
df = pd.read_csv(INPUT)
required = {"species", "reaction", "metabolite", "coefficient"}

if not required <= set(df.columns):
    raise ValueError(f"Input must contain columns: {sorted(required)}")

if PROJECTION not in ("shared_metabolite", "substrate_product"):
    raise ValueError("Unknown PROJECTION rule.")

if NULL_MODEL not in ("bipartite", "projection"):
    raise ValueError("Unknown NULL_MODEL.")

groups = df.groupby("species")["group"].first() if "group" in df.columns else None
df = df.dropna(subset=list(required))
df["coefficient"] = pd.to_numeric(df["coefficient"], errors="coerce")
df = df.dropna(subset=["coefficient"])
df = df.groupby(["species", "reaction", "metabolite"], as_index=False)["coefficient"].sum()
df = df[df["coefficient"] != 0]

if df.empty:
    raise ValueError("No valid reaction-metabolite records.")


def project(inc, reactions):
    """Reaction-node projection of the reaction-metabolite incidence table."""
    G = nx.Graph()
    G.add_nodes_from(reactions)

    for _, g in inc.groupby("metabolite"):
        if PROJECTION == "shared_metabolite":
            G.add_edges_from(combinations(g["reaction"].unique(), 2))
        else:
            producers = g.loc[g["coefficient"] > 0, "reaction"].unique()
            consumers = g.loc[g["coefficient"] < 0, "reaction"].unique()
            G.add_edges_from((p, c) for p in producers for c in consumers if p != c)

    return G


def topology(G):
    """L, C and r on the largest connected component (unreachable pairs never enter L)."""
    lcc = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    n = lcc.number_of_nodes()

    if n < 3 or lcc.number_of_edges() == 0:
        return lcc, np.nan, np.nan, np.nan

    d = shortest_path(nx.to_scipy_sparse_array(lcc), unweighted=True, directed=False)
    L = d.sum() / (n * (n - 1))

    # Nodes with degree < 2 have an undefined clustering coefficient and are excluded
    local = [c for node, c in nx.clustering(lcc).items() if lcc.degree(node) >= 2]
    C = np.mean(local) if local else np.nan

    return lcc, L, C, nx.degree_assortativity_coefficient(lcc)


def bipartite_null(inc):
    """Degree-preserving swaps of reaction-metabolite incidences (keeps each reaction's roles)."""
    r = inc["reaction"].to_numpy()
    m = inc["metabolite"].to_numpy().copy()
    pairs = set(zip(r, m))

    for _ in range(SWAPS_PER_EDGE * len(r)):
        i, j = rng.integers(len(r), size=2)

        if r[i] == r[j] or m[i] == m[j]:
            continue
        if (r[i], m[j]) in pairs or (r[j], m[i]) in pairs:
            continue

        pairs -= {(r[i], m[i]), (r[j], m[j])}
        pairs |= {(r[i], m[j]), (r[j], m[i])}
        m[i], m[j] = m[j], m[i]

    return pd.DataFrame({"reaction": r, "metabolite": m, "coefficient": inc["coefficient"].to_numpy()})


# 2. Analyse each species, with and without currency-metabolite removal
settings = {"all metabolites": []}
if CURRENCY:
    settings["currency removed"] = CURRENCY

rows = []

for species, sp in df.groupby("species"):
    reactions = sp["reaction"].unique()

    for setting, drop in settings.items():
        inc = sp[~sp["metabolite"].isin(drop)]
        G = project(inc, reactions)
        lcc, L, C, r = topology(G)
        sizes = sorted((len(c) for c in nx.connected_components(G)), reverse=True)

        # 3. Randomized reference networks
        Lr, Cr, rr = [], [], []
        if lcc.number_of_nodes() >= 4 and G.number_of_edges() >= 4:
            for _ in range(N_RANDOM):
                if NULL_MODEL == "bipartite":
                    Gr = project(bipartite_null(inc), reactions)
                else:
                    Gr = G.copy()
                    nswap = SWAPS_PER_EDGE * G.number_of_edges()
                    nx.double_edge_swap(Gr, nswap=nswap, max_tries=100 * nswap,
                                        seed=int(rng.integers(1_000_000_000)))
                _, l, c, a = topology(Gr)
                Lr.append(l)
                Cr.append(c)
                rr.append(a)

        Lr, Cr, rr = np.array(Lr, float), np.array(Cr, float), np.array(rr, float)

        with np.errstate(all="ignore"):
            sigma = (C / np.nanmean(Cr)) / (L / np.nanmean(Lr)) if Lr.size else np.nan
            per_run = (C / Cr) / (L / Lr) if Lr.size else np.array([np.nan])

        rows.append({
            "species": species,
            "group": groups.get(species) if groups is not None else np.nan,
            "setting": setting,
            "reactions": len(reactions),
            "metabolites": inc["metabolite"].nunique(),
            "edges": G.number_of_edges(),
            "density": nx.density(G),
            "components": len(sizes),
            "lcc_nodes": lcc.number_of_nodes(),
            "lcc_fraction": sizes[0] / len(reactions),
            "efficiency": nx.global_efficiency(G),
            "L": L, "C": C, "r": r,
            "L_random": np.nanmean(Lr) if Lr.size else np.nan,
            "C_random": np.nanmean(Cr) if Cr.size else np.nan,
            "r_random": np.nanmean(rr) if rr.size else np.nan,
            "sigma": sigma,
            "sigma_2.5%": np.nanpercentile(per_run, 2.5),
            "sigma_97.5%": np.nanpercentile(per_run, 97.5),
        })

res = pd.DataFrame(rows)
res.to_csv("topology_by_species.csv", index=False)

print(f"Projection: {PROJECTION} | null model: {NULL_MODEL} ({N_RANDOM} networks) | seed {SEED}")
print("L, C and r are computed on the largest connected component; efficiency uses all reactions.")
print(res.drop(columns=["group"]).round(4).to_string(index=False))

# 4. Cross-species comparison (all-metabolite networks)
base = res[res["setting"] == "all metabolites"].dropna(subset=["L", "C", "r"]).copy()

if len(base) >= 4:
    print("\nSpearman correlations across species (association only, no phylogenetic control):")
    for a, b in [("lcc_nodes", "L"), ("L", "C"), ("L", "r")]:
        rho, p = stats.spearmanr(base[a], base[b])
        print(f"  {a} vs {b}: rho = {rho:.3f}, p = {p:.3f}")

if len(base) >= 3:
    # Remove the dependence of L on network size before comparing groups
    fit = stats.linregress(np.log(base["lcc_nodes"]), base["L"])
    base["L_size_adjusted"] = base["L"] - (fit.intercept + fit.slope * np.log(base["lcc_nodes"]))
    print(f"\nL vs ln(LCC size): slope = {fit.slope:.3f}, p = {fit.pvalue:.3f}")
    print("Coefficient of variation of L across species:", round(base["L"].std() / base["L"].mean(), 3))

    samples = [g["L_size_adjusted"].to_numpy() for _, g in base.groupby("group")] if groups is not None else []
    if len(samples) >= 2:
        try:
            H, p = stats.kruskal(*samples)
            print(f"Kruskal-Wallis, size-adjusted L across groups: H = {H:.3f}, p = {p:.3f}")
        except ValueError as err:
            print("Group comparison not possible:", err)

# 5. Figure: observed vs randomized L, C and r per species
fig, axes = plt.subplots(1, 3, figsize=(13, 4))
for ax, (obs, rnd, label) in zip(axes, [("L", "L_random", "average path length L"),
                                       ("C", "C_random", "average clustering C"),
                                       ("r", "r_random", "degree assortativity r")]):
    ax.bar(base["species"], base[obs], label="observed")
    ax.plot(base["species"], base[rnd], "k_", markersize=18, label="randomized mean")
    ax.set_ylabel(label)
    ax.tick_params(axis="x", rotation=60)
axes[0].legend()
plt.tight_layout()
plt.savefig("small_world_assortativity.png", dpi=200)
