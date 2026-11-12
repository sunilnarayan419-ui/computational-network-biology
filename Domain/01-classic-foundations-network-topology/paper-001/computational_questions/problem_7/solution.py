import numpy as np
import pandas as pd
import networkx as nx
import powerlaw
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Configuration
INPUT = "edges.csv"          # columns: source,target
DIRECTED = False             # True -> in-degree and out-degree analysed separately
SIGNIFICANCE = 0.10          # bootstrap goodness-of-fit threshold (reject power law if p < this)
LR_ALPHA = 0.05              # significance for likelihood-ratio comparisons
BOOTSTRAPS = 200             # semi-parametric bootstrap replicates
MIN_NODES = 50               # minimum network size for a classification
MIN_TAIL = 50                # minimum number of observations in the fitted tail
ALTERNATIVES = ["exponential", "lognormal", "truncated_power_law"]
SEED = 42

np.random.seed(SEED)         # powerlaw draws from NumPy's global generator
rng = np.random.default_rng(SEED)

# 1. Load edge list and build graph (self-loops removed, duplicate edges merged)
edges = pd.read_csv(INPUT, usecols=["source", "target"]).dropna().astype(str)
edges = edges[edges["source"] != edges["target"]]

if edges.empty:
    raise ValueError("No valid edges after removing missing values and self-loops.")

G = nx.from_pandas_edgelist(
    edges, "source", "target",
    create_using=nx.DiGraph if DIRECTED else nx.Graph
)

print("Nodes:", G.number_of_nodes(), "| Edges:", G.number_of_edges())
print("Software: powerlaw", powerlaw.__version__, "| numpy", np.__version__, "| seed", SEED)

# 2. Degree sequences (degree-zero nodes cannot enter a discrete power-law tail)
if DIRECTED:
    measures = {"in-degree": dict(G.in_degree()), "out-degree": dict(G.out_degree())}
else:
    measures = {"degree": dict(G.degree())}

fig, axes = plt.subplots(len(measures), 2, figsize=(10, 4 * len(measures)), squeeze=False)

for row, (name, degrees) in enumerate(measures.items()):
    k = np.array(list(degrees.values()))
    k = k[k > 0]

    table = pd.Series(k).value_counts().sort_index().rename_axis("degree").reset_index(name="nodes")
    table["P(k)"] = table["nodes"] / k.size
    table.to_csv(f"degree_frequency_{name}.csv", index=False)

    print(f"\n=== {name} ===")
    if np.unique(k).size < 3:
        print("Too few distinct degree values for distribution fitting: Inconclusive")
        continue

    # 3. Maximum-likelihood power-law fit with KS-based k_min selection
    fit = powerlaw.Fit(k, discrete=True, verbose=False)
    alpha, kmin, ks_stat = fit.alpha, fit.xmin, fit.D
    n_tail = int((k >= kmin).sum())

    # 4. Semi-parametric bootstrap goodness of fit (Clauset et al., 2009)
    tail, body = k[k >= kmin], k[k < kmin]
    p_tail = tail.size / k.size
    exceed, valid = 0, 0

    for _ in range(BOOTSTRAPS):
        n_from_tail = rng.binomial(k.size, p_tail)
        synthetic = [fit.power_law.generate_random(n_from_tail)] if n_from_tail else []

        if body.size and k.size - n_from_tail:
            synthetic.append(rng.choice(body, k.size - n_from_tail))

        try:
            refit = powerlaw.Fit(np.concatenate(synthetic).astype(int), discrete=True, verbose=False)
            exceed += refit.D >= ks_stat
            valid += 1
        except Exception:
            continue

    p_gof = exceed / valid if valid else np.nan

    # 5. Likelihood-ratio comparison against alternatives (R > 0 favours the power law)
    comparisons = []
    for alt in ALTERNATIVES:
        R, p = fit.distribution_compare(
            "power_law", alt, nested=(alt == "truncated_power_law"), normalized_ratio=True
        )
        comparisons.append((alt, R, p, getattr(fit, alt).parameters))

    # 6. Evidence-based classification
    favoured = [a for a, R, p, _ in comparisons if R < 0 and p < LR_ALPHA]
    beats_exponential = any(a == "exponential" and R > 0 and p < LR_ALPHA for a, R, p, _ in comparisons)

    if G.number_of_nodes() < MIN_NODES or n_tail < MIN_TAIL:
        verdict = "Inconclusive (network or fitted tail too small)"
    elif not np.isfinite(p_gof):
        verdict = "Inconclusive (bootstrap failed)"
    elif p_gof < SIGNIFICANCE:
        verdict = "Non-scale-free evidence (power law rejected by goodness of fit)"
    elif favoured:
        verdict = "Non-scale-free evidence (alternative preferred: " + ", ".join(favoured) + ")"
    elif "exponential" in ALTERNATIVES and not beats_exponential:
        verdict = "Inconclusive (light-tailed exponential model not excluded)"
    else:
        verdict = "Scale-free supported (power law not rejected)"

    # 7. Report
    print(f"Nodes with degree > 0: {k.size}")
    print(f"alpha = {alpha:.3f} +/- {fit.sigma:.3f} | k_min = {kmin:g} | tail size = {n_tail}")
    print(f"KS statistic D = {ks_stat:.4f} | bootstrap p = {p_gof:.3f} ({valid} valid replicates)")

    for alt, R, p, params in comparisons:
        print(f"Power law vs {alt}: normalized R = {R:.3f}, p = {p:.4f}, parameters = {params}")

    print("Classification:", verdict)
    print("Note: failing to reject a power law is not proof that the network is scale-free.")

    # 8. Degree distribution (log-log) and CCDF
    ax = axes[row, 0]
    ax.loglog(table["degree"], table["P(k)"], "o")
    ax.set(xlabel="degree k", ylabel="P(k)", title=f"{name} distribution")

    ax = axes[row, 1]
    fit.plot_ccdf(ax=ax, original_data=True, label="data")
    fit.power_law.plot_ccdf(ax=ax, linestyle="--", label="power law")
    for alt in ALTERNATIVES:
        getattr(fit, alt).plot_ccdf(ax=ax, linestyle=":", label=alt)
    ax.set(xlabel="degree k", ylabel="CCDF", title=f"{name} CCDF (k_min = {kmin:g})")
    ax.legend(fontsize=7)

plt.tight_layout()
plt.savefig("degree_distribution.png", dpi=200)
