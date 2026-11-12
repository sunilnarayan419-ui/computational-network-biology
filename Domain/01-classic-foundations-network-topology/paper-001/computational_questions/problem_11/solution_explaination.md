# Solution 11: Hierarchical Modularity in Cellular Metabolic Networks

## 1. Solution Overview

This solution looks for hierarchical modular organization in a cellular metabolic network built from a stoichiometric model. It finds modules at several resolutions, measures how well they are separated, how well fine modules sit inside coarser ones, how stable they are, and how they relate to curated metabolic subsystems. Each signature is compared with randomized metabolic networks.

The analysis is restricted to cellular metabolic networks supplied as SBML models. The reaction–metabolite relationships are the primary data and are kept as a bipartite incidence table. The reaction and metabolite networks used for module detection are projections of that table.

The implementation follows six computational stages:

1. Reading the model and building the reaction–metabolite incidence.
2. Reaction and metabolite projections.
3. Multi-resolution module detection and signature metrics.
4. Module connectivity and parent–child relationships.
5. Randomized metabolic networks.
6. Signature summary, sensitivity, and figure.

## 2. Technology Selection

| Library | Computational responsibility |
|---|---|
| COBRApy | Reading the SBML model, reaction stoichiometry, subsystem annotation, boundary-reaction identification |
| Pandas | Incidence table, cross-tabulations, result tables |
| NumPy | Randomization and summary statistics |
| NetworkX | Weighted graphs, Louvain community detection, modularity |
| scikit-learn | Normalized mutual information (NMI) between partitions |
| Matplotlib | Modularity and nestedness figure |

NetworkX's built-in Louvain implementation is used, so `igraph` and `leidenalg` are not needed.

## 3. Computational Methodology

### Stage 1: Incidence from the Stoichiometric Model

**Library:** COBRApy, Pandas

Every non-boundary reaction contributes one row per metabolite with its stoichiometric coefficient (negative for substrates, positive for products). Exchange, demand, and sink reactions are excluded because they connect only one metabolite. Biomass or other reactions can be removed with `EXCLUDE_REACTIONS`. Metabolite identifiers keep their compartment suffix, so compartment-specific metabolites stay distinct. The filtered incidence is saved as `bipartite_incidence.csv`.

Reaction subsystems come from the model's own annotation. Reactions without a subsystem are left unannotated and are not assigned one.

### Stage 2: Network Representations

**Library:** NetworkX

Two projections are analysed.

| Representation | Nodes | Edge between two nodes when | Weight |
|---|---|---|---|
| Reaction | Reactions | They share a metabolite (any role) | Number of shared metabolites |
| Metabolite | Metabolites | One reaction converts one into the other (substrate to product) | Number of reactions doing so |

Currency metabolites can create large artificial hubs. `CURRENCY` is empty by default, so nothing is removed unless the user supplies identifiers. When it is supplied, the observed analysis is also repeated with all metabolites retained so the effect can be compared. Isolated nodes are not part of the projected graphs.

### Stage 3: Multi-Resolution Modules and Signature Metrics

**Library:** NetworkX, scikit-learn

Modules are detected with weighted Louvain community detection at each resolution $\gamma$ in `RESOLUTIONS`. Larger $\gamma$ gives smaller modules. The partition from the first run is used for descriptions, and `N_RUNS` runs with different seeds are used to measure stability.

Modularity of a partition, calculated at the standard resolution so values are comparable:

$$Q = \frac{1}{2m}\sum_{ij}\left(A_{ij} - \frac{k_i k_j}{2m}\right)\delta(c_i, c_j)$$

where $A_{ij}$ is the (weighted) adjacency, $k_i$ the node strength, $m$ the total edge weight, and $\delta$ equals 1 when nodes $i$ and $j$ are in the same module.

Four quantities are computed at each resolution:

| Quantity | Definition |
|---|---|
| Modularity $Q$ | As above |
| Stability | Mean NMI between the partitions of repeated runs |
| Nestedness | For consecutive resolutions, the fraction of nodes that lie in the coarse module that best overlaps their fine module: $\sum_f \max_c \lvert f \cap c\rvert / N$ |
| Subsystem agreement | NMI between modules and curated subsystem labels (reaction network, annotated reactions only) |

Nestedness is measured and not forced. Louvain partitions at different resolutions need not be perfectly nested, and curated subsystems are not assumed to form a perfect hierarchy. Reactions that share metabolites across subsystems are never split off or reassigned.

### Stage 4: Connectivity and Parent–Child Tables

**Library:** Pandas

For every module at every resolution, the script reports its size, the weight of internal edges, the weight of edges to other modules, and the internal fraction:

$$f_{\text{int}} = \frac{w_{\text{internal}}}{w_{\text{internal}} + w_{\text{external}}}$$

For consecutive resolutions, each fine module is assigned to the coarse module holding most of its nodes, with the fraction that falls in that parent. Fractions below 1 show overlap and are not hidden. Files written: `module_connectivity.csv`, `module_parent_child.csv`, and `module_assignments_<representation>.csv` (module labels at every resolution plus subsystem).

### Stage 5: Randomized Metabolic Networks

**Library:** NumPy, Pandas

`N_RANDOM` randomized networks are generated by swapping pairs of reaction–metabolite incidences. Each reaction keeps its number of substrates and products, and each metabolite keeps its number of reactions. Randomization is applied to the bipartite table before projection, so the projections of randomized networks have the same structural origin as the observed ones. The complete pipeline is rerun on each randomized network.

For each metric, the script reports the observed value, the null mean and standard deviation, a z-score, and the empirical one-sided p-value:

$$p = \frac{1 + \#\{\text{null} \ge \text{observed}\}}{1 + N_{\text{random}}}$$

### Stage 6: Signature Summary and Figure

For each representation the script lists which of the following screening signatures are met:

1. Modularity above the null at every resolution ($p \le 0.05$).
2. Nestedness at least `NESTED_MIN` at every pair of consecutive resolutions and above the null mean.
3. Stability at least `STABILITY_MIN` at every resolution.
4. Subsystem agreement above the null at one or more resolutions (reaction network only).

The thresholds are configurable screening values, not established constants. A positive summary is evidence consistent with hierarchical modularity, not proof, and no single signature is sufficient. The figure `hierarchical_modularity.png` plots observed and randomized $Q$ and nestedness against resolution, and the full table is written to `hierarchical_modularity_report.csv`.

## 4. Input Data Requirements

| Item | Requirement |
|---|---|
| `MODEL` | SBML file of a cellular metabolic reconstruction with reaction stoichiometry |
| Subsystem annotation | Optional but needed for the subsystem comparison |
| `CURRENCY`, `EXCLUDE_REACTIONS` | Optional user-supplied identifier lists |

Only cellular metabolic networks are supported. Protein-interaction and gene-regulatory networks are out of scope.

## 5. Expected Computational Output

| Output | Description |
|---|---|
| Model summary | Reactions, metabolites, compartments, annotation coverage |
| Report table | Observed and null values for $Q$, stability, nestedness, subsystem NMI |
| Signature list | Which screening signatures are met per representation |
| Module tables | Connectivity, parent–child, and assignment files |
| Currency sensitivity | Observed metrics with and without the listed metabolites |
| `hierarchical_modularity.png` | Figure of $Q$ and nestedness |

## 6. Scientific Assumptions and Limitations

1. Results depend on the projection rule, currency handling, and the quality and completeness of the reconstruction.
2. Louvain is heuristic and resolution-dependent, and different resolutions are not guaranteed to form a hierarchy. Nestedness is therefore reported as a measured quantity.
3. Modularity can be high in randomized networks with heterogeneous degrees, which is why the null comparison and several signatures are required.
4. Subsystem labels are curated by different groups and are not ground truth, and unannotated reactions are excluded from that comparison.
5. The bipartite null preserves degrees and substrate–product role counts but not other biochemical constraints such as mass balance or reaction chemistry.
6. The script does not test whether each hierarchical level is scale-free, and does not assume it should be. Scale-free analysis of module-level degree distributions is not implemented.
7. Hierarchical-clustering and dendrogram methods, directed projections, and pathway-level enrichment tests are not implemented.
8. The metabolite projection ignores reaction reversibility and counts substrate–product pairs as written.

## 7. Conclusion

The script builds reaction and metabolite networks from a stoichiometric model while retaining the underlying reaction–metabolite incidence. It detects modules at several resolutions, measures their quality, stability, nesting, and subsystem agreement, and compares each with randomized metabolic networks. Hierarchical modularity is judged from several signatures together and never from one resolution alone.
