# Solution 13: Vulnerability Assessment and Damage-Minimizing Resilience of Pyruvate Reaction Pipelines Under Unbuffered or Weakly Buffered Conditions

## 1. Solution Overview

This solution simulates a pyruvate-to-product reaction pipeline in an unbuffered or weakly buffered medium, applies defined disturbances, measures the resulting loss of performance, and evaluates user-defined mitigation strategies.

pH is not obtained by adding up produced protons. It is solved from electroneutrality (charge balance) of the supplied medium, including weak acids with their pKa values, strong ions, and activity corrections, and it is recalculated at every time step as the reaction changes the composition. Pathway performance is computed from a stoichiometric model with COBRApy, with enzyme capacities set by measured response curves.

All biochemical inputs are supplied by the user. The script contains no built-in medium composition, buffer capacity, pKa, or enzyme parameter.

The implementation follows seven computational stages:

1. Loading and validating inputs.
2. Stoichiometric model construction.
3. Baseline medium and performance.
4. Per-reaction vulnerability screen.
5. Disturbance scenarios and mitigation strategies.
6. Comparison of baseline, disturbed, and mitigated performance.
7. Figure.

## 2. Technology Selection

| Library | Computational responsibility |
|---|---|
| COBRApy | Stoichiometric model and flux optimisation at each time step |
| Pandas | Input tables, results, and summaries |
| NumPy | Speciation arithmetic and time courses |
| SciPy | Root finding for pH |
| Matplotlib | Damage and pH figure |

NetworkX, SALib, eQuilibrator, and a dedicated chemical-equilibrium package are not used. The medium chemistry needed here is a charge balance of listed weak acids and strong ions, which is solved directly. Sensitivity is explored through scenario magnitude ranges and not through a SALib method.

## 3. Computational Methodology

### Stage 1: Input Data

**Library:** Pandas

| Table | Contents |
|---|---|
| Stoichiometry | Reaction, metabolite, coefficient |
| Reactions | Reference capacity `vmax_M_per_h` (mol/L/h) and optional reversibility |
| Response | Measured relative activity of each enzyme as a function of pH, temperature ($T_K$), or ionic strength ($I_M$) |
| Medium | Species, kind (acid, neutral, strong cation, strong anion), concentration, charge, pKa values (separated by `;` for polyprotic acids) |
| Scenarios | Disturbance kind, magnitude, optional low and high magnitude, time window, target reaction, weight |
| Mitigations | Kind, value, optional target, cost |

For acids, charge is that of the fully protonated form. For strong ions, charge is the absolute charge. Species such as the substrate and product must appear in the medium table with their acid properties. Unlisted species are not part of the acid–base calculation. Validation stops the run if the medium is missing or has an unknown species kind, so pH is never guessed.

### Stage 2: Stoichiometric Model

**Library:** COBRApy

A COBRApy model is built from the stoichiometry. Boundary reactions allow supply of the source (limited by what is left in the medium), removal of the target, and free exchange of any metabolites in `EXCHANGE`. The objective is to maximise the export flux of the target. The network must be mass-balanced for internal species, so cofactors must be regenerated in the supplied reactions or listed in `EXCHANGE`.

### Stage 3: Medium Chemistry and Time-Stepped Simulation

**Library:** SciPy, NumPy, COBRApy

For each time step, the script does the following.

**pH from electroneutrality.** With strong-ion difference $\text{SID} = \sum z_{\text{cation}} c - \sum z_{\text{anion}} c$ plus any added strong base, the pH satisfies:

$$\text{SID} + [\mathrm{H^+}] - [\mathrm{OH^-}] + \sum_i C_i\, \bar q_i(\mathrm{pH}) = 0$$

where $C_i$ is the total concentration of weak acid $i$ across all protonation states and $\bar q_i$ is its mean charge. For an acid with fully protonated charge $z_0$ and successive pKa values, the ratio of consecutive states is:

$$\frac{[A_j]}{[A_{j-1}]} = 10^{\,\mathrm{pH} - \mathrm{pK}_{a,j}} \exp\left[-A_\gamma\,\bigl(2(z_0 - j) + 1\bigr)\right]$$

with $[\mathrm{H^+}] = 10^{-\mathrm{pH}} e^{A_\gamma}$ and $[\mathrm{OH^-}] = 10^{\mathrm{pH} - \mathrm{pK}_w} e^{A_\gamma}$. The activity term comes from the extended Debye–Hückel expression:

$$A_\gamma = \frac{\alpha(T)\sqrt{I}}{1 + B\sqrt{I}}, \quad \alpha(T) = 1.10708 - 1.54508\times10^{-3}\,T + 5.95584\times10^{-6}\,T^2, \quad B = 1.6$$

The ionic strength $I = \tfrac{1}{2}\sum c_i z_i^2$ depends on speciation, so pH and $I$ are iterated to convergence. The equation has one root between pH 0 and 14, found with Brent's method. If no root exists the run stops with an error and reports that the composition is insufficient. The same function gives the baseline pH and a numerical buffer capacity, $\beta = \Delta(\text{base added}) / \Delta\text{pH}$, calculated by adding and removing a small amount of strong base.

Because pH comes from the total concentrations of the listed acid–base species, a reaction that converts one acid into another changes the pH through that equilibrium, and no proton counter is used.

**Capacities.** Each reaction's capacity is its reference `vmax` multiplied by the measured relative activity at the current pH, temperature, and ionic strength, interpolated linearly inside the measured range. Factors for different variables are multiplied, which assumes they act independently. A response value outside the measured range uses the edge value and is counted and reported. A reaction with no response data is not assumed insensitive. It is listed as unassessed, and it enters the numbers as unaffected.

**Flux and composition.** The model is optimised under these capacities, and the resulting fluxes update the totals of the medium species over the time step. The substrate supply is limited by the amount left.

**Performance** is the amount of target formed over the horizon:

$$J = C_{\text{target}}(t_{\text{end}}) - C_{\text{target}}(0)$$

**Damage** is the relative loss against the undisturbed baseline:

$$D = 1 - \frac{J_{\text{scenario}}}{J_{\text{baseline}}}$$

Negative $D$ means the disturbance increased performance, since not every pH change is harmful.

### Stage 4: Vulnerability Screen

**Library:** COBRApy

Each reaction's activity is scaled by `VULN_FACTOR` over the whole horizon, one reaction at a time, and $D$ is recorded. Reactions with the highest $D$ are the pathway bottlenecks under this perturbation. The table is saved as `vulnerability_screen.csv` and records whether response data were available.

### Stage 5: Disturbance Scenarios and Mitigation

| Scenario kind | Effect |
|---|---|
| `acid_base_pulse` | Strong base (positive) or acid (negative) added from the start time, in mol/L |
| `ionic_salt` | Inert 1:1 salt added from the start time, which raises ionic strength without changing charge balance |
| `T_shift` | Temperature change in kelvin during the window |
| `substrate_scale` | Factor on available substrate during the window |
| `enzyme_scale` | Factor on one enzyme's activity (or all) during the window |
| `data_corrupt` | Response data of a reaction (or all) scaled by the magnitude, representing erroneous input |

Biochemical scenarios are labelled `biochemical`. The `data_corrupt` scenario is labelled `data_integrity` because it represents faulty input to the pipeline and not a change in the chemistry. Exceptions raised by the software are recorded as software failures, separately from performance loss.

| Mitigation kind | Effect |
|---|---|
| `ph_control` | Holds pH at a set value |
| `T_control` | Holds temperature at a set value |
| `enzyme_boost` | Multiplies enzyme activity |
| `add_buffer` | Adds a buffer species already defined in the medium table |
| `data_validation` | Rejects response data that are negative or above `REL_MAX` |

Every scenario is run at its low, nominal, and high magnitude with every single mitigation and every combination up to `MAX_COMBO`. Mitigations are only evaluated. The script never decides which are chemically or operationally allowed, so users must supply only realistic options with their costs. Data validation can only catch values outside its range, and corrupted values that stay in range are not detected.

### Stage 6: Comparison and Selection

For each scenario the output gives baseline, unmitigated, and mitigated damage with its low–high range, plus:

- **Residual damage:** loss in production rate averaged over the last tenth of the horizon.
- **Recovery time:** time after the end of the disturbance window (or after the start, for persistent disturbances) until the production rate returns to within `RECOVERY_TOL` of baseline. It is empty if recovery does not happen within the horizon.

The best mitigation per scenario minimises:

$$\text{objective} = D_{\text{nominal}} + \lambda \cdot \text{cost}$$

where $\lambda$ is `COST_WEIGHT`. Across all scenarios, strategies are ranked by weighted expected damage plus cost, and the cheapest strategy whose worst-case damage over all scenarios and magnitudes is at most `MAX_DAMAGE` is reported.

### Stage 7: Figure

**Library:** Matplotlib

`pipeline_resilience.png` shows damage per scenario with and without the best mitigation, and the pH trajectory of the most damaging scenario compared with the baseline and the mitigated case.

## 4. Input Data Requirements

See Stage 1. Rate units must be mol/L/h and concentrations mol/L. The pKa values and $\mathrm{p}K_w$ are assumed valid at `PKA_T_K`, and the script warns when the baseline temperature is more than 5 K away. No value is estimated or downloaded by the script.

## 5. Expected Computational Output

| Output | Description |
|---|---|
| Baseline pH, ionic strength, buffer capacity, $J_0$ | Medium and performance baseline |
| `vulnerability_screen.csv` | Damage per reaction when its activity is scaled |
| `scenario_mitigation_results.csv` | Damage (low, nominal, high), residual damage, recovery, cost for each scenario and strategy |
| `strategy_summary.csv` | Expected and worst-case damage per strategy |
| Console summary | Baseline vs disturbed vs best mitigated, best strategies, software failures, out-of-range response use |
| `pipeline_resilience.png` | Figure |

## 6. Scientific Assumptions and Limitations

1. pH is only as good as the medium specification. Missing buffering species (for example dissolved CO$_2$ or unlisted cofactors) make the pH prediction incomplete, and the script can warn only when a measured initial pH is supplied.
2. Pyruvate and the other species follow the supplied pKa values with extended Debye–Hückel activity corrections, which are approximate at higher ionic strength. Metal binding and temperature-dependent pKa values are not modelled.
3. Enzyme response to pH, temperature, and ionic strength comes only from the supplied curves, and the effects are assumed independent and multiplicative.
4. A pH controller is modelled as holding the pH exactly at its setpoint, and the ionic strength is still taken from the uncontrolled calculation.
5. Flux comes from optimisation of a steady-state model at each step, so regulation and kinetics beyond the supplied capacities are not represented.
6. Thermodynamic feasibility is not re-evaluated here, and Problem 12 deals with it.
7. Scenario ranges give only a simple bracket of uncertainty, and no global sensitivity analysis is performed.
8. Runtime grows with the number of scenarios, mitigation combinations, and time steps, because each requires one optimisation per step.
9. Results are computational predictions and are not experimental proof.

## 7. Conclusion

The script combines an equilibrium description of an unbuffered medium with a stoichiometric pathway model, defines performance loss explicitly, separates biochemical vulnerability from data-integrity and software failures, and evaluates user-defined mitigations for damage, cost, residual loss, and recovery. All conclusions depend on the supplied medium and enzyme-response data.
