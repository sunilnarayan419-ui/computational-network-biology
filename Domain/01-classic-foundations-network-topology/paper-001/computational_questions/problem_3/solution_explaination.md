# Solution 03: Scale-Free Metabolic Network and Pyruvate Hub Identification

## 1. Solution Overview

This solution builds a metabolite network from user-supplied biochemical reactions, calculates the connectivity of every metabolite, and reports pyruvate's actual degree, rank and centrality. It then examines the degree distribution with a statistical power-law fit and compares the network with degree-preserving randomized networks.

Pyruvate is not assumed to be a hub. Its connectivity is calculated from the data, and the power-law analysis reports evidence and its limits rather than a "scale-free" verdict.

The script follows seven stages:

1. Reaction data loading.
2. Metabolite network construction.
3. Connectivity and centrality calculation.
4. Pyruvate analysis.
5. Degree distribution and power-law analysis.
6. Degree-preserving randomization.
7. Visualization.

## 2. Technology Selection

| Library | Computational responsibility |
|---|---|
| Pandas | Reading reactions, metabolite table, ranking, CSV output |
| NumPy | Numerical summaries and empirical p-values |
| NetworkX | Directed and undirected graphs, centralities, degree-preserving edge swaps |
| powerlaw | Maximum-likelihood power-law fit and likelihood-ratio comparison with alternative distributions |
| Matplotlib | Degree-distribution plot and hub sub-network figure |

SciPy is not imported directly. The `powerlaw` package handles the fitting statistics.

## 3. Computational Methodology

### Stage 1: Reaction Data Loading

**Library:** Pandas

The input has one row per reaction with `reaction_id`, `substrates`, `products` and an optional `reversible` column. Multiple metabolites in one cell are separated by `;`. Duplicate reaction IDs are removed and counted. If `reversible` is absent, all reactions are treated as irreversible, and this assumption is stated in the code.

### Stage 2: Metabolite Network Construction

**Library:** NetworkX

Multi-substrate and multi-product reactions are represented as follows:

- Every substrate is connected by a directed edge to every product of the same reaction (substrate → product).
- For reversible reactions, the opposite edges are added as well.
- Substrates of the same reaction are not connected to each other, and products are not connected to each other.

Currency metabolites (`CURRENCY`, default water, proton, ATP, ADP, phosphate, NAD(H), NADP(H)) are removed before the edges are created when `EXCLUDE_CURRENCY` is true. Compartment suffixes such as `_c` are ignored when matching names. Coenzyme A is not on the default list, so it is retained and can be examined as a hub. The list can be changed in the configuration.

The directed graph `D` is the primary representation. An undirected copy `U` (a link exists if either direction exists) is used for degree, degree centrality, betweenness and randomization. Stoichiometric coefficients are not used, because only the topology is analysed.

### Stage 3: Connectivity and Centrality

**Library:** Pandas, NetworkX

The table contains, for each metabolite:

| Column | Definition |
|---|---|
| `degree` | k(i) = Σ A(i,j): number of distinct neighbouring metabolites in `U` |
| `in_degree` / `out_degree` | Distinct metabolites converted into / from this metabolite in `D` |
| `n_reactions` | Number of distinct reactions in which it appears |
| `produced_by` / `consumed_by` | Number of reactions producing / consuming it; reversible reactions count in both |
| `degree_centrality` | k(i) / (N − 1) |
| `betweenness` | Fraction of shortest paths between other metabolite pairs passing through it (undirected graph; exact unless `BETWEENNESS_K` is set) |

Reaction participation and distinct-neighbour degree are kept as separate columns because they are not equivalent. The table is sorted by degree and saved as `metabolite_metrics.csv`; the top `N_HUBS` are printed.

### Stage 4: Pyruvate Analysis

**Library:** Pandas

Identifiers in `PYRUVATE_IDS` are matched case-insensitively against the data. Different compartments (for example `pyr_c` and `pyr_m`) are different nodes and are reported separately if listed. For each match the script prints its degree, rank, the percentage of metabolites with a smaller degree, all table values, the reaction IDs and the neighbouring metabolites. If no identifier is found, the script says so.

### Stage 5: Degree Distribution and Power-Law Analysis

**Library:** powerlaw, Pandas

The degree distribution is

P(k) = N(k) / N

where N is the number of metabolites, and it is plotted on log-log axes for total, in- and out-degree (only k > 0 can be shown).

For the power-law model P(k) ∝ k^(−γ), the script uses `powerlaw.Fit` with `discrete=True`. The package estimates the lower cut-off x_min by minimizing the Kolmogorov–Smirnov distance and estimates γ by maximum likelihood for the degrees ≥ x_min. The script prints γ with its standard error, x_min, the KS distance and the number of nodes in the tail.

The power law is compared with exponential, log-normal and truncated power-law distributions by a likelihood-ratio test. R > 0 favours the power law, and the p-value indicates whether the sign of R is reliable.

The script does not run a goodness-of-fit bootstrap and does not output a "scale-free" decision. A linear-looking log-log plot is not used as evidence.

### Stage 6: Degree-Preserving Randomization

**Library:** NetworkX, NumPy

The undirected network is rewired `N_RANDOM` times with `double_edge_swap`, using `SWAPS_PER_EDGE × (number of edges)` swaps per realization. Each swap exchanges the end points of two edges, so every metabolite keeps its degree. The average clustering coefficient and degree assortativity of the observed network are compared with the randomized values. The empirical p-value is

p = (1 + number of random values with |x − mean| ≥ |observed − mean|) / (N_RANDOM + 1)

This test is skipped for networks with fewer than 4 nodes or 10 edges. Because degrees are preserved, the degree distribution itself is not compared; Problem 04 handles the Erdős–Rényi comparison.

### Stage 7: Visualization

**Library:** Matplotlib, NetworkX

One PNG contains the log-log degree distributions on the left and, on the right, the sub-network induced by the `N_HUBS` highest-degree metabolites with pyruvate highlighted in gold and node size scaled with degree. It is a hub sub-network, not the full network.

## 4. Input Data Requirements

| reaction_id | substrates | products | reversible |
|---|---|---|---|
| R001 | glucose | pyruvate | false |
| R002 | pyruvate;coa;nad | acetyl_coa;nadh;co2 | false |

- The table should preferably come from a curated organism-specific reconstruction.
- Metabolite identifiers must be consistent, and compartment-specific identifiers are treated as separate nodes.
- `PYRUVATE_IDS` must use the identifiers in the file.

The rows above only illustrate the format. No dataset is bundled, and no result is provided.

## 5. Output Interpretation

| Output | Interpretation |
|---|---|
| Hub table | Metabolites with the most distinct neighbours; high degree alone does not prove biological importance |
| Pyruvate rank and percentage | Where pyruvate lies in this dataset; it is a hub only if the numbers show a very high rank |
| Betweenness | How often a metabolite lies on shortest connecting paths in the network |
| alpha, x_min, tail size | Degree exponent and the part of the distribution it describes; a small tail gives an unreliable estimate |
| R and p against alternatives | Relative support for a power law versus other heavy-tailed or light-tailed models |
| Randomization p-values | Whether clustering or assortativity differ from networks with the same degrees |

## 6. Assumptions and Limitations

1. Network topology depends strongly on the reconstruction, the identifier mapping and the handling of currency metabolites.
2. The substrate → product representation is one of several possible projections. Shortest paths and betweenness follow this definition and ignore stoichiometry, enzyme constraints and thermodynamics.
3. Undirected betweenness and degree-preserving randomization ignore reaction direction.
4. Metabolic networks of limited size give a short tail, so the power-law fit and the model comparison have low power. Without a goodness-of-fit test, a non-rejected power law is not evidence that the network is scale-free.
5. A high degree does not by itself show hub function in cells; pathway activity and flux are not evaluated.
6. The script does not compare pyruvate's reactions with an external database; this validation must be done separately.
7. Compartmentalized identifiers split a metabolite into several nodes, which lowers each node's degree.
