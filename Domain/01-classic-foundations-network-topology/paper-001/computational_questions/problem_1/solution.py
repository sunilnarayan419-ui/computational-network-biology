import numpy as np
import pandas as pd
import networkx as nx

from scipy.stats import t
from sklearn.impute import KNNImputer
from statsmodels.stats.multitest import multipletests

# Configuration
INPUT = "microarray.csv"
SOURCE = "GENE_A"
TARGETS = ["GENE_B", "GENE_C"]
FDR_LIMIT = 0.05
MIN_CORRELATION = 0.5

# 1. Load and preprocess expression matrix
df = pd.read_csv(INPUT, index_col=0)
X = df.apply(pd.to_numeric, errors="coerce").to_numpy()

if X.shape[1] < 3:
    raise ValueError("At least 3 samples are required.")

if not np.isfinite(X).all():
    X = KNNImputer(n_neighbors=5).fit_transform(X.T).T

if not np.isfinite(X).all():
    raise ValueError("Unresolved missing values.")

# 2. Calculate pairwise correlations and p-values
R = np.corrcoef(X)
n = X.shape[1]

upper = np.triu_indices_from(R, k=1)
r = R[upper]

valid = np.isfinite(r) & (np.abs(r) > 0)
r = r[valid]

# Pearson correlation significance test
dfree = n - 2
stat = np.abs(r) * np.sqrt(dfree / (1 - r**2))
p = 2 * t.sf(stat, dfree)

# 3. Multiple-testing correction
q = multipletests(p, method="fdr_bh")[1]

i, j = upper
i, j = i[valid], j[valid]

# 4. Construct weighted cellular network
G = nx.Graph()
G.add_nodes_from(df.index)

for a, b, corr, fdr in zip(i, j, r, q):
    if fdr <= FDR_LIMIT and abs(corr) >= MIN_CORRELATION:
        G.add_edge(
            df.index[a],
            df.index[b],
            weight=-np.log(abs(corr)),
            correlation=corr,
            fdr=fdr
        )

# 5. Identify shortest paths
for target in TARGETS:
    try:
        path = nx.shortest_path(
            G, SOURCE, target, weight="weight"
        )
        cost = nx.shortest_path_length(
            G, SOURCE, target, weight="weight"
        )

        probability_score = np.exp(-cost)

        print(f"\n{SOURCE} -> {target}")
        print("Path:", " -> ".join(path))
        print("Path cost:", round(cost, 4))
        print("Combined association score:", round(probability_score, 6))

    except (nx.NetworkXNoPath, nx.NodeNotFound):
        print(f"\nNo valid path: {SOURCE} -> {target}")

# 6. Network summary
print("\nNetwork nodes:", G.number_of_nodes())
print("Network edges:", G.number_of_edges())
