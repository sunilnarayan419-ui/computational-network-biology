# Solution 05: Dynamic Scale-Free Network Modelling with Biochemical Feasibility

## 1. Solution Overview

This solution is a single-machine prototype that searches for pyruvate-associated pathways in a stoichiometric metabolic reconstruction while evaluating the feasibility of every reaction direction. Pathways are found by a cost-weighted k-shortest-path search on a reaction–metabolite graph, then checked for stoichiometric consistency with a linear program and, as a reference, compared with flux balance analysis (FBA).

A topological path is never treated as proof of feasibility. Missing biochemical evidence is recorded as missing and given a penalty, whereas violated constraints exclude a reaction direction.

The script follows six stages:

1. Loading the reaction store.
2. Reaction feasibility evaluation.
3. Dynamic updates with version history and a preferential-attachment diagnostic.
4. Feasibility assessment of the final network.
5. Pathway search, stoichiometric check and FBA reference.
6. Unconstrained baseline.

## 2. Technology Selection

| Library | Computational responsibility |
|---|---|
| COBRApy | Reading the SBML model; stoichiometry, bounds, gene rules; FBA and parsimonious FBA reference |
| Pandas | Thermodynamic and concentration tables, feasibility table, ranked output |
| NumPy | Stoichiometric matrix, logarithms |
| NetworkX | Bipartite search graph and k-shortest simple paths |
| SciPy | Linear program for stoichiometric consistency of a candidate pathway |

The eQuilibrator API is **not** called by the script. Transformed standard Gibbs energies depend on pH, ionic strength, pMg and temperature, so the user supplies ΔG°′ values (for example, calculated with eQuilibrator under chosen conditions) in `thermo.csv`. No Gibbs energies or enzyme parameters are created by the script. Distributed infrastructure is not implemented.

## 3. Computational Methodology

### Stage 1: Reaction Store

**Library:** COBRApy

The SBML model is read and every non-boundary reaction is stored with its full stoichiometry (`metabolite → coefficient`), lower and upper flux bounds, and an enzyme status. The enzyme status is evaluated from the gene–reaction rule: genes in `UNAVAILABLE_GENES` are set to false and `and` / `or` are evaluated logically. It is `True` (enzyme available), `False` (required enzyme unavailable) or `None` (no gene annotation). Boundary (exchange, demand, sink) reactions are not pathway steps.

Optional tables:

- `thermo.csv`: `reaction_id`, `dG0_prime_kJ_mol` and optionally `dG0_sd_kJ_mol`, for the reaction as written in the model.
- `concentrations.csv`: `metabolite`, `conc_min_M`, `conc_max_M`.

### Stage 2: Reaction Feasibility F(r)

**Library:** NumPy, Pandas

Each reaction is assessed separately in the forward (`fwd`) and reverse (`rev`) direction.

**Directionality (D).** A direction is blocked when its flux bounds do not allow flux (upper bound ≤ 0 for forward, lower bound ≥ 0 for reverse).

**Substrate and cofactor availability (S, C).** A direction is rejected if it consumes a metabolite in `UNAVAILABLE_METABOLITES`. Metabolites in `POOL` (currency and cofactors) are not exchanged between reactions in the search but are subject to this rule.

**Enzyme (E).** A direction is rejected when the gene rule evaluates to false. If no annotation exists, the evidence is recorded as missing. Kinetic parameters and enzyme abundances are not used.

**Thermodynamics (T).**

ΔG′ = ΔG°′ + RT ln Q, with Q = Π c_i^(ν_i)

(ν_i > 0 for products, < 0 for substrates, R = 8.314 × 10⁻³ kJ mol⁻¹ K⁻¹, T = `TEMPERATURE_K`). Because concentrations are ranges, an interval is calculated:

- ΔG°′ is widened by `Z_SD` × its standard deviation if one is given.
- For each product, ln c contributes its smallest value at c_min and its largest at c_max; for each substrate the roles are reversed.
- Protons and water are left out of Q, because ΔG°′ is a transformed value.
- Metabolites without an entry in `concentrations.csv` use `CONC_DEFAULT` (1 µM–10 mM). This is an explicit assumption and not measured data.
- Reactions that involve more than one compartment (transport) are treated as having no thermodynamic data, because the membrane potential and pH gradient are not supplied.

For the reverse direction, the interval becomes [−ΔG′_max, −ΔG′_min]. If the whole interval is above zero, the direction is rejected as thermodynamically unfavourable. If the whole interval is below zero, it is favourable. If it contains zero, the direction is retained with a weight equal to the fraction of the interval below zero (at least `W_MIN`). If ΔG°′ is missing, the direction is retained and the missing evidence is recorded.

**Score and cost.** Each missing evidence type (enzyme annotation, thermodynamic data) multiplies F by `W_UNKNOWN`; a partly uncertain thermodynamic interval multiplies F by the weight described above. The cost is

C(r) = −ln F(r),  with 0 < F(r) ≤ 1

These are ordinal evidence-completeness weights and not calibrated probabilities. They should be calibrated against curated evidence before quantitative interpretation. Rejected directions have F = 0 and infinite cost.

### Stage 3: Dynamic Updates and Preferential-Attachment Diagnostic

**Library:** NumPy, standard library

If `UPDATES_CSV` is given, rows are applied in order. `add` inserts a reaction with its stoichiometry and bounds (enzyme status unknown; its ΔG°′ can be supplied in `thermo.csv`), `remove` deletes one and `set_bounds` changes its bounds. After every update, an independent copy of the whole store is saved together with a fingerprint (SHA-256 of reaction IDs, stoichiometry and bounds). At the end, all stored versions are re-fingerprinted to confirm that earlier versions were not modified.

For every added reaction, the script prints a constrained preferential-attachment diagnostic for the metabolites it uses:

Π(i) = k(i) B(i) / Σ_j k(j) B(j)

k(i) is the number of reactions containing metabolite i before the update, and B(i) is the biochemical compatibility factor, defined here as 1 for non-pool metabolites and 0 for pool metabolites. Π(i) shows how much attachment probability the growth hypothesis assigns to each metabolite. This is a network-growth hypothesis and not a proven evolutionary mechanism. The script does not generate synthetic networks and does not test the hypothesis statistically. Heterogeneous connectivity is reported as the most connected metabolites (reaction participation, pool excluded).

### Stage 4: Feasibility of the Final Network

The final reaction store is evaluated and written to `reaction_feasibility.csv` (admissibility, F, cost, ΔG′ interval, rejection reasons, missing evidence). The screen output lists how many directions are admissible, how many are rejected with a reason (directions blocked only by flux bounds are not listed) and how many admissible directions have missing evidence. Missing evidence is therefore distinguished from infeasibility.

### Stage 5: Pathway Search, Stoichiometric Check and FBA Reference

**Library:** NetworkX, SciPy, COBRApy

**Search graph.** A directed bipartite graph is built from admissible reaction directions. Metabolite nodes connect to reaction nodes (cost C(r) plus a small `TIE_BREAK`), and reaction nodes connect to their product nodes (cost 0). Pool metabolites are not nodes. Reaction identity is kept, so alternative reactions between the same metabolites are different paths. `shortest_simple_paths` returns up to `K_PATHS` simple (non-repeating) paths from `SOURCE` to `TARGET` in order of cost; paths with more than `MAX_LEN` reactions are discarded. The ranked cost reported for a pathway is ΣC(r) without the tie-break.

**Stoichiometric check.** A multi-substrate reaction needs all its substrates, which a graph path alone does not ensure. For each candidate, a linear program tests whether the oriented reactions can run at steady state while producing one unit of `TARGET`:

- Flux v_j ≥ 0 for each pathway reaction (in its direction), limited by its bound and `FLUX_CAP`.
- For every metabolite that is not in the pool and is not `SOURCE`, net production (S v) ≥ 0, except `TARGET`, for which S v ≥ 1.
- The objective minimizes the total flux.

If the program is infeasible, the script lists co-substrates that the pathway consumes but never produces. Pool metabolites and `SOURCE` are assumed to be supplied, so cofactor regeneration is not tested here.

**Ranking.** Pathways that pass the check are ranked first, then by cost and length (`ranked_pathways.csv`). The file also lists missing evidence, ΔG′ intervals per reaction and unmet co-substrates.

**FBA reference.** On the unmodified input model, a demand reaction for `TARGET` is added and maximized, optionally with the medium in `FBA_MEDIUM`. If feasible, parsimonious FBA gives a flux distribution, and the script reports how many pathway reactions carry flux in the pathway direction. This comparison is not applied to dynamic updates and depends on the medium used. FBA may produce `TARGET` from other sources than pyruvate.

### Stage 6: Unconstrained Baseline

The shortest path with unit cost per reaction, using only flux-bound directionality, is found and printed. For each of its reactions, the script shows whether it is admissible under the feasibility rules and why not, and whether the path passes the stoichiometric check. This shows how a purely topological path can differ from the feasibility-aware result.

## 4. Input Data Requirements

| Input | Format |
|---|---|
| Model | SBML file with metabolite IDs (compartment suffix `_c`, `_m` ...), bounds and gene–reaction rules |
| Thermodynamics (optional) | `reaction_id`, `dG0_prime_kJ_mol`, optional `dG0_sd_kJ_mol` (kJ/mol) |
| Concentrations (optional) | `metabolite`, `conc_min_M`, `conc_max_M` (mol/L) |
| Updates (optional) | `action`, `reaction_id`, `stoichiometry` as `met:coef;met:coef`, `lower_bound`, `upper_bound` |

Metabolite identifiers in the configuration (`SOURCE`, `TARGET`, `POOL`, unavailable sets) must match the model's identifiers. No data are bundled with the repository.

## 5. Output Interpretation

| Output | Interpretation |
|---|---|
| Network versions and digests | History of updates; the intact-history check confirms no earlier version was altered |
| `reaction_feasibility.csv` | Evidence-based assessment per reaction direction, with reasons and missing evidence |
| `ranked_pathways.csv` | Candidate pathways with cost, stoichiometric support, missing evidence and FBA activity |
| Low cost | Little missing or uncertain evidence along the path, not a higher flux or probability |
| Unmet co-substrates | Substrates that the path needs but does not make |
| FBA activity | Agreement with a constraint-based model under its medium |
| Unconstrained baseline | Illustration of the topological shortest path and its feasibility problems |

## 6. Assumptions and Limitations

1. Thermodynamic favourability does not establish flux. Enzyme kinetics and abundance are not modelled.
2. Concentration ranges that are missing are replaced by the default range, which affects ΔG′ intervals. Results should be read together with `missing_evidence`.
3. Feasibility weights are ordinal and not calibrated probabilities.
4. Pool metabolites are assumed available and regenerated. Redox and energy balance of the whole cell are not checked by the path analysis; FBA gives only a partial control.
5. The path search follows main substrate → product relations; the linear program is applied afterwards. Searches are limited to `K_PATHS` candidates and `MAX_LEN` reactions.
6. Transport reactions have no thermodynamic evaluation because membrane potential and pH gradients are not supplied.
7. The preferential-attachment diagnostic only evaluates attachment propensities of added reactions and does not support any claim about metabolic network evolution.
8. The scripts do not plot networks or run distributed computation. The prototype is a single process.
9. FBA is run on the unmodified input model; dynamic updates are not applied to it.
