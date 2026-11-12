# Solution 08: Small-World Effect and Degree Assortativity in Evolving Reaction Networks

## 1. Solution Overview

This solution measures average path length, clustering, small-world character, and degree assortativity in reaction-node networks built from metabolic reaction data, and compares these properties across species or evolutionary groups.

The reaction–metabolite incidence table is the primary data representation. A reaction-node network is derived from it using an explicit projection rule. Small-world character is assessed against randomized reference networks, and the same metrics are then compared across species while controlling for network size.

The implementation follows six computational stages:

1. Loading and validating reaction–metabolite data.
2. Reaction-node projection.
3. Topology metrics on the largest connected component.
4. Randomized reference networks and the small-world index.
5. Currency-metabolite sensitivity run.
6. Cross-species comparison and visualization.

## 2. Technology Selection

| Library | Computational responsibility |
|---|---|
| Pandas | Loading incidence data and collecting results |
| NumPy | Random numbers and summary statistics |
| NetworkX | Graph construction, clustering, assortativity, components, global efficiency |
| SciPy | All-pairs shortest paths (`csgraph`) and statistical tests |
| Matplotlib | Comparison figure |

COBRApy is not used. The script reads the reaction–metabolite table directly, so users who start from a genome-scale model must export stoichiometry to that table first.

## 3. Computational Methodology

### Stage 1: Data Loading

**Library:** Pandas

Each row of the input states that a reaction in a given species involves a metabolite with a stoichiometric coefficient. Negative coefficients are substrates and positive coefficients are products. Rows with missing values or non-numeric coefficients are dropped, repeated records are summed, and records whose net coefficient is zero are removed. An optional `group` column assigns each species to an evolutionary group.

Metabolite identifiers should already encode compartment (for example `atp_c`), so compartment-specific metabolites remain distinct.

### Stage 2: Reaction-Node Projection

**Library:** NetworkX

Reactions are nodes. The edge rule is set by `PROJECTION`:

| Rule | Edge between reactions A and B when |
|---|---|
| `shared_metabolite` | A and B involve at least one common metabolite, in any role |
| `substrate_product` | A produces a metabolite that B consumes, or the reverse |

The resulting graph is undirected and unweighted. The substrate–product rule uses the coefficient signs as written and does not treat reversible reactions as running backwards. Reactions that share no metabolite with any other reaction remain as isolated nodes.

### Stage 3: Topology Metrics

**Library:** SciPy, NetworkX

Unreachable reaction pairs have infinite distance, so average path length is calculated only on the largest connected component (LCC). The script also reports the number of components, the LCC fraction, and the global efficiency of the whole network, which handles unreachable pairs without discarding them.

Average shortest-path length on the LCC with $N$ reaction nodes:

$$L = \frac{1}{N(N-1)} \sum_{i \ne j} d(i,j)$$

Local clustering of node $i$ with degree $k_i$ and $e_i$ edges among its neighbours:

$$C_i = \frac{2 e_i}{k_i (k_i - 1)}$$

Nodes with $k_i < 2$ have an undefined $C_i$ and are excluded from the mean $C$.

Degree assortativity is the Pearson correlation between the degrees at the two ends of an edge:

$$r = \operatorname{corr}(k_{\text{source}}, k_{\text{target}})$$

Positive $r$ means high-degree reactions tend to link to high-degree reactions, negative $r$ means they link to low-degree reactions. This is degree assortativity only. It says nothing about functional or pathway similarity.

### Stage 4: Randomized Reference and Small-World Index

**Library:** NumPy, NetworkX

Two reference models are available through `NULL_MODEL`:

| Model | Procedure |
|---|---|
| `bipartite` (default) | Reaction–metabolite incidences are swapped in pairs so every reaction keeps its number of substrates and products and every metabolite keeps its number of reactions. The reaction projection is then recalculated. |
| `projection` | Degree-preserving edge swaps are applied directly to the projected reaction network. |

The bipartite null is the default because a projection contains many triangles by construction, since every metabolite produces a clique of reactions. Rewiring the projection destroys these triangles and would inflate the small-world index for trivial reasons.

The small-world index is:

$$\sigma = \frac{C / C_{\text{random}}}{L / L_{\text{random}}}$$

where $C_{\text{random}}$ and $L_{\text{random}}$ are means over `N_RANDOM` randomized networks. The script also reports the 2.5th and 97.5th percentiles of $\sigma$ calculated per random network, as a measure of uncertainty. Randomized metrics are also taken on each randomized network's LCC.

$\sigma > 1$ is a descriptive indicator and does not establish small-world structure alone.

### Stage 5: Currency-Metabolite Sensitivity

If `CURRENCY` lists metabolite identifiers, the full analysis is repeated with those metabolites removed before projection. Both versions appear in the output table so the influence of currency metabolites can be inspected. The script never removes any metabolite unless the user supplies the list.

### Stage 6: Cross-Species Comparison

**Library:** SciPy, Matplotlib

For networks with all metabolites retained, the script reports:

- Spearman correlations between LCC size and $L$, between $L$ and $C$, and between $L$ and $r$.
- A linear fit of $L$ on $\ln(N_{\text{LCC}})$, whose residuals give a size-adjusted $L$.
- The coefficient of variation of $L$ across species, as a simple measure of stability.
- A Kruskal–Wallis test of size-adjusted $L$ across evolutionary groups, when a `group` column with at least two groups is present.

The figure `small_world_assortativity.png` shows observed values with the randomized means for $L$, $C$ and $r$. The full table is saved as `topology_by_species.csv`.

## 4. Input Data Requirements

| Column | Meaning |
|---|---|
| `species` | Species or evolutionary-state identifier |
| `reaction` | Reaction identifier |
| `metabolite` | Metabolite identifier, ideally with compartment |
| `coefficient` | Stoichiometric coefficient (negative substrate, positive product) |
| `group` (optional) | Evolutionary group of the species |

Species must come from comparable reconstructions. Different annotation standards and database versions are not corrected by the script.

## 5. Expected Computational Output

| Output | Description |
|---|---|
| Node, edge, density, component statistics | Network size and connectivity per species |
| $L$, $C$, $r$ | Observed topology metrics on the LCC |
| Randomized means and $\sigma$ with percentile range | Small-world comparison |
| Global efficiency | Path-based measure covering disconnected networks |
| Cross-species statistics | Correlations, size-adjusted $L$, group test |
| `topology_by_species.csv`, `small_world_assortativity.png` | Saved results |

Interpretation: if $L$ stays similar across species while network size changes, the data are consistent with a stable path length. If $L$ tracks size, the change is explained by size and not necessarily by anything evolutionary.

## 6. Scientific Assumptions and Limitations

1. Metrics depend on the projection rule, and the shared-metabolite rule can give very different networks from the substrate–product rule.
2. Currency metabolites can create artificial connectivity, and their handling is left to the user and tested only through the optional sensitivity run.
3. Species are treated as independent, but related species are phylogenetically non-independent, so the tests are exploratory.
4. Differences in reconstruction quality and completeness are not modelled, and size adjustment is only a rough control.
5. Topology alone cannot show that evolution maintains or selects a particular path length.
6. Edge removal for incomplete data, alternative projection comparisons, and phylogenetic comparative methods are not implemented.
7. Path lengths are computed with a dense distance matrix, so memory grows with the square of the LCC size.

## 7. Conclusion

The script quantifies path length, clustering, small-world character, and degree assortativity in reaction networks derived from reaction–metabolite data, and compares them across species with size adjustment and randomized references. The results describe network organization and do not by themselves establish an evolutionary mechanism.
