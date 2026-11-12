# Solution 06: Hierarchical Dynamic Metabolic Network Modelling

## 1. Solution Overview

This solution extends the Problem 05 concept with hierarchical metabolic organization. Reactions and their stoichiometry remain the primary biochemical representation. Four levels are derived from them and from a reaction-membership table: subsystems (Level 1), pathways (Level 2), reactions (Level 3) and metabolites (Level 4).

Pathways are explored from coarse to fine. Routes between subsystems select candidate subsystems, routes between pathways inside them select candidate pathways, and the reaction-level search is then limited to the reactions of the selected pathways. The same search is also run without the hierarchy, and the two results are compared.

Reaction feasibility is reused from Problem 05 through its output table instead of being recalculated. The script does not assume that any level is scale-free.

The script follows six stages:

1. Loading reactions and hierarchy memberships.
2. Reusing reaction costs from Problem 05.
3. Dynamic updates with consistency records.
4. Construction and description of the four levels.
5. Coarse-to-fine and non-hierarchical pathway search.
6. Visualization.

## 2. Technology Selection

| Library | Computational responsibility |
|---|---|
| COBRApy | Reading the SBML model (stoichiometry, bounds, subsystem annotation) |
| Pandas | Hierarchy memberships, feasibility table, results |
| NumPy | Degree statistics, stoichiometric matrix |
| NetworkX | Bipartite reaction–metabolite graph, projections, group graphs, k-shortest paths |
| SciPy | Linear program for stoichiometric consistency |
| powerlaw | Power-law fit and comparison with a log-normal model, where enough nodes exist |
| Matplotlib | Level 1 and Level 2 network figure |

FBA is not repeated here; it is available in Problem 05. Distributed computation is not implemented.

## 3. Computational Methodology

### Stage 1: Reactions and Hierarchy

**Library:** COBRApy, Pandas

Every non-boundary reaction is stored with its stoichiometry and flux bounds. The hierarchy is read from `hierarchy.csv` (`reaction_id`, `subsystem`, `pathway`). A reaction can appear in several rows, so overlapping pathway membership is kept. A pathway is identified by `subsystem | pathway`, so the same pathway name in different subsystems stays separate.

If no hierarchy file is given, the SBML subsystem is used for both Level 1 and Level 2 (pathway = subsystem). The script states this; it does not invent pathways.

Before cleaning, the script reports memberships that refer to unknown reactions, duplicated membership rows and reactions without annotation. Duplicates and unknown reactions are then removed. Unannotated reactions remain in the model but are not reachable through the hierarchy.

### Stage 2: Reaction Costs

**Library:** Pandas

For each reaction direction that the flux bounds allow, the cost is taken from `reaction_feasibility.csv`, produced by the Problem 05 script:

C(r) = −ln F(r)

- Directions listed as not admissible are excluded.
- Directions with no row in the table (for example reactions added by an update) receive `UNKNOWN_COST`. Missing evidence is therefore handled as uncertainty and not as infeasibility.
- If no feasibility table is given, all bound-allowed directions have `UNKNOWN_COST`, and only flux-bound directionality constrains the search.

The table must be calculated for the same model and conditions.

### Stage 3: Dynamic Updates

**Library:** Pandas, NumPy

Updates in `UPDATES_CSV` are applied sequentially: `add` (stoichiometry, bounds, subsystem and pathway required) or `remove`. The hierarchy is derived from the reaction store and the membership table, so a change reaches all levels when the levels are built. After the initial state and each update, a record lists the number of reactions, subsystems and pathways and the consistency counts (unknown reactions, duplicate memberships, unannotated reactions).

For an added reaction, the script also prints a hierarchical preferential-attachment diagnostic:

Π(i) ∝ k(i) B(i) M(i)

k(i) is the number of reactions containing metabolite i, B(i) = 1 for non-pool metabolites and 0 for pool metabolites (pool metabolites are not included), and M(i) = 1 if metabolite i occurs in a reaction of the added reaction's own pathway, otherwise `M_OTHER`. Π is normalized over all existing metabolites. `M_OTHER` is a hypothesis parameter, not an estimated value. The diagnostic is not tested against data and says nothing about biological evolution.

### Stage 4: Four Levels

**Library:** NetworkX, NumPy, powerlaw

| Level | Object | Edges |
|---|---|---|
| 4 / 3 backbone | Bipartite graph of reactions and non-pool metabolites with the coefficient on each edge | Reaction contains metabolite |
| 3 | Reaction network (bipartite projection) | Two reactions share a non-pool metabolite |
| 2 | Pathway network | Two pathways have reactions that share a non-pool metabolite |
| 1 | Subsystem network | Two subsystems have reactions that share a non-pool metabolite |

Cross-level relationships are explicit: a reaction belongs to a pathway and subsystem through the membership table, and a group is connected to another only through metabolites that the reactions of both actually share. No other shortcuts are created. Shared metabolites and overlapping memberships are preserved. Highly connected non-pool metabolites can link many groups, so the pool list should be reviewed for the dataset.

For each level, the script prints the number of nodes, the mean and maximum degree and ⟨k²⟩/⟨k⟩ (a heterogeneity measure). Only when a level has at least `MIN_FIT_NODES` nodes and at least five distinct degrees, a discrete power-law fit is made with `powerlaw` (exponent γ, x_min and a likelihood-ratio comparison with a log-normal). Otherwise, the script states that no scale-free claim can be made. Each level is assessed independently. Goodness-of-fit bootstrap and randomized reference networks are not run here; Problems 03 and 04 provide those tools.

### Stage 5: Pathway Search

**Library:** NetworkX, SciPy

**Search graph and check.** As in Problem 05, a bipartite directed graph of reaction directions and non-pool metabolites is built with edge costs C(r) plus `TIE_BREAK`, and `K_PATHS` lowest-cost simple paths are retrieved from `SOURCE` to `TARGET`. Each candidate is checked by a linear program with the same assumptions as in Problem 05 (pool metabolites and `SOURCE` supplied, intermediates not net consumed, one unit of `TARGET` produced). Unmet co-substrates are reported. The bottleneck reaction of a pathway is its reaction with the highest cost, and the subsystems and pathways traversed are listed.

**Coarse-to-fine exploration.**

1. *Subsystems:* start subsystems contain reactions with `SOURCE`, and end subsystems contain reactions with `TARGET`. Up to `K_COARSE` shortest routes in the subsystem network are found for each start/end pair, and all subsystems on these routes are selected.
2. *Pathways:* the pathway network restricted to the selected subsystems is searched in the same way, and the pathways on the routes are selected.
3. *Reactions:* the reaction-level search is limited to the reactions that belong to the selected pathways.

Coarse routes use hop counts only, because feasibility costs are defined for reactions. The hierarchical result is therefore a heuristic that can miss the best pathway.

**Non-hierarchical search.** The same search runs on all reactions. The script compares the number of edges in the search graph, the run time, the number of candidates, the best cost among stoichiometrically supported pathways and the cost gap between the two approaches. Both result sets are written to `hierarchical_pathways.csv`, with a `search` column. The most connected metabolites (reactions per metabolite) are printed as a hub summary.

### Stage 6: Visualization

**Library:** Matplotlib

One figure shows the subsystem network and the pathway network of the selected subsystems, with selected nodes in orange. Reaction and metabolite levels are described by statistics and not plotted.

## 4. Mathematical Operations

**Cost of a pathway:** C(π) = Σ C(r_i), with C(r) = −ln F(r).

**Heterogeneity:** ⟨k²⟩/⟨k⟩, where k is degree. It increases with a broad degree distribution and is not by itself evidence for a power law.

**Power-law exponent:** P(k) ∝ k^(−γ), estimated by maximum likelihood for degrees ≥ x_min by the `powerlaw` package.

**Stoichiometric condition:** S v ≥ 0 for intermediates, (S v)_TARGET ≥ 1, 0 ≤ v ≤ bound, solved as a linear program.

## 5. Input Data Requirements

| Input | Format |
|---|---|
| Model | SBML file, metabolite IDs with compartment suffix |
| Feasibility table | Output of Problem 05: `reaction_id`, `direction`, `admissible`, `cost` |
| Hierarchy | `reaction_id`, `subsystem`, `pathway`; several rows per reaction are allowed |
| Updates (optional) | `action`, `reaction_id`, `stoichiometry` (`met:coef;met:coef`), `lower_bound`, `upper_bound`, `subsystem`, `pathway` |

`SOURCE`, `TARGET` and `POOL` must use the model's identifiers. No data are bundled with the repository.

## 6. Output Interpretation

| Output | Interpretation |
|---|---|
| Hierarchy consistency and update log | Missing, duplicated or inconsistent relationships between levels |
| Level statistics | Connectivity at each level; the power-law fit is reported only when enough nodes exist |
| Coarse search summary | How much of the reaction set the hierarchy keeps |
| Comparison table | Cost, size and time of hierarchical and non-hierarchical searches |
| Cost gap | A positive value means the hierarchical search missed a cheaper supported pathway |
| `hierarchical_pathways.csv` | Candidate pathways with cost, stoichiometric support, bottleneck, subsystems and pathways |

## 7. Assumptions and Limitations

1. A low cost reflects few missing or uncertain feasibility items, not a higher flux or probability. Topological connectivity, hierarchical membership, thermodynamic feasibility and flux are different properties.
2. The hierarchy depends on the annotation file. With only SBML subsystems, Level 2 repeats Level 1.
3. The coarse stage uses topology only and can exclude a better pathway. The reported gap measures this for the given data.
4. Run-time comparisons are only informative for networks large enough that search size matters.
5. The feasibility table must match the model, conditions and updates used here. Updated reactions without a table row get `UNKNOWN_COST`.
6. Pool metabolites are assumed available and regenerated, and cofactor balance of the whole cell is not tested. FBA is not repeated.
7. Power-law results for the small subsystem and pathway levels are not meaningful, and the script does not make them.
8. The preferential-attachment diagnostic with `B` and `M` is a hypothesis and its parameters are not tested against data.
9. Compartment-specific metabolites are separate nodes, and compartment compatibility is checked only through the identifiers in each reaction.
