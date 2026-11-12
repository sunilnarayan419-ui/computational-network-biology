# Solution 12: Robustness of Pyruvate Reaction Pathways Under Variable Temperature and Salt Conditions

## 1. Solution Overview

This solution asks whether pathways from pyruvate to a chosen product remain thermodynamically feasible when temperature and ionic conditions change. For each condition it recalculates the Gibbs free energy of every pathway reaction, converts the ionic composition into ionic strength, applies activity and temperature corrections only where the supplied data allow, and compares reaction feasibility and pathway ranking across conditions.

All biochemical values (standard energies, enthalpies, charges, concentrations, ionic composition) are supplied by the user. The script contains no built-in reaction energies, and missing information is reported as missing evidence and never as infeasibility.

The implementation follows six computational stages:

1. Loading and validating user data.
2. Pathway enumeration.
3. Ionic strength and condition-dependent $\Delta G$.
4. Reaction robustness and sensitivity.
5. Pathway ranking across conditions.
6. Visualization.

## 2. Technology Selection

| Library | Computational responsibility |
|---|---|
| Pandas | Input tables and result tables |
| NumPy | Thermodynamic arithmetic and Monte Carlo draws |
| SciPy | Normal distribution, Kendall rank correlation |
| NetworkX | Enumerating pathways in a metabolite–reaction graph |
| Matplotlib | Heat map of $\Delta G$ by condition |

COBRApy, eQuilibrator, and SALib are not used. The script evaluates thermodynamics only, so no flux model is needed. Values from eQuilibrator or the literature can be supplied through the input tables, and no sensitivity-analysis method from SALib is implemented.

## 3. Computational Methodology

### Stage 1: Input Data

**Library:** Pandas

| Table | Contents |
|---|---|
| Reactions | $\Delta_r G'^{\circ}$ (kJ/mol) with its reference temperature, ionic strength, and pH; optional standard deviation and $\Delta_r H'^{\circ}$ |
| Stoichiometry | Reaction, metabolite, coefficient (negative substrate, positive product) |
| Metabolites | Charge and concentration (mol/L) |
| Conditions | Temperature (K) and pH; the first row is the baseline |
| Ions | Ion concentration and charge for each condition |

Each reaction's $\Delta_r G'^{\circ}$ is valid only at its stated reference conditions. Species whose activity is fixed at 1 in the transformed energy (such as water and protons) are listed in `UNIT_ACTIVITY` and left out of the reaction quotient and activity corrections. Temperatures must be in kelvin.

### Stage 2: Pathway Enumeration

**Library:** NetworkX

A graph links metabolites to reaction steps, where each reaction can be traversed in either direction. Currency metabolites listed in `CURRENCY` are not used as links, which prevents trivial routes through cofactors, but they still enter the energy calculation. All simple paths from `SOURCE` to `TARGET` up to `MAX_PATH_LEN` reactions are collected, up to `MAX_PATHS`. A path that uses the same reaction twice is discarded. Thermodynamics decides whether each step is favorable in the direction traversed.

### Stage 3: Ionic Strength and $\Delta G$

**Library:** NumPy

Salt concentration and ionic strength are different quantities. Ionic strength is calculated from the supplied ionic composition:

$$I = \frac{1}{2}\sum_i c_i z_i^2$$

where $c_i$ is the molar concentration and $z_i$ the charge of ion $i$. Only supplied ions contribute. A condition without ion rows has $I = 0$.

Activity coefficients use the extended Debye–Hückel form with the temperature-dependent parameter of Alberty:

$$\ln \gamma_i = -\frac{\alpha(T)\, z_i^2 \sqrt{I}}{1 + B\sqrt{I}}, \qquad \alpha(T) = 1.10708 - 1.54508\times10^{-3}\,T + 5.95584\times10^{-6}\,T^2, \quad B = 1.6$$

These are physical-chemistry model constants fixed in the script. They are not biochemical measurements, and they are an approximation valid at low to moderate ionic strength.

The condition-dependent reaction energy is built in three steps, with absolute temperature $T$ throughout.

1. Remove the effect of the reference ionic strength:

$$\Delta G^{\circ}_{I=0} = \Delta G'^{\circ}_{\text{ref}} - R\,T_{\text{ref}} \sum_i \nu_i \ln \gamma_i(I_{\text{ref}}, T_{\text{ref}})$$

2. If $\Delta H'^{\circ}$ is supplied, move to temperature $T$ with the Gibbs–Helmholtz relation, assuming a constant enthalpy of reaction:

$$\Delta G^{\circ}(T) = \Delta G^{\circ}_{I=0}\,\frac{T}{T_{\text{ref}}} + \Delta H^{\circ}\left(1 - \frac{T}{T_{\text{ref}}}\right)$$

3. Add the activity-based reaction quotient at condition $(T, I)$:

$$\Delta G = \Delta G^{\circ}(T) + R\,T \sum_i \nu_i \left[\ln \gamma_i(I, T) + \ln c_i\right]$$

where $\nu_i$ are the stoichiometric coefficients (positive for products) and $c_i$ the supplied concentrations in mol/L.

Data gaps are handled explicitly:

| Situation | Result |
|---|---|
| Missing $\Delta G'^{\circ}$ or concentration | $\Delta G$ is not calculated and the reaction is marked missing evidence |
| Condition pH differs from the reaction's reference pH (beyond `PH_TOL`) | Not calculated, since no pH-correction data are supplied |
| Missing charge | Calculated assuming an ideal solution, with a note |
| Missing $\Delta H'^{\circ}$ | Only the $RT \ln Q$ term responds to temperature, with a note |

Each condition is processed in its own error-handling block, so a failing condition is recorded in the status table and does not stop the other conditions.

### Stage 4: Reaction Feasibility and Sensitivity

**Library:** SciPy, NumPy

If a standard deviation of $\Delta G'^{\circ}$ is supplied, the probability that a reaction runs forward is $P(\Delta G < 0) = \Phi(-\Delta G / \sigma)$. Without a standard deviation, the sign of $\Delta G$ is used. A reaction is labelled forward if $P \ge$ `P_THRESHOLD`, reverse if $P \le 1 - $ `P_THRESHOLD`, and otherwise near equilibrium or uncertain.

At the baseline condition, one-factor-at-a-time sensitivities are obtained as the change in $\Delta G$ for a temperature change of $\pm$`SENS_DT` and an ionic-strength change of $\pm$`SENS_DI`. A reaction is flagged condition-sensitive if either change reaches `SENSITIVE_KJ`, or if it is forward in one condition and reverse in another.

### Stage 5: Pathway Ranking

**Library:** NumPy, SciPy

For a pathway traversed in a given direction, the step energy is $\pm \Delta G$. The pathway is described by its **bottleneck**, the least favorable step:

$$\Delta G_{\text{bottleneck}} = \max_{s}\, \Delta G_s$$

A more negative bottleneck means a larger minimum driving force. The total energy is also reported. The probability that all steps are feasible at once is estimated by Monte Carlo sampling of $\Delta G'^{\circ}$ errors, using the same random numbers for every condition. Errors are drawn independently for different reactions, which is an assumption. Pathways with any missing-evidence step are labelled undetermined. Kendall rank correlation compares the pathway ranking at the baseline with each other condition.

### Stage 6: Visualization

**Library:** Matplotlib

The figure `reaction_dG_by_condition.png` is a heat map of $\Delta G$ for each pathway reaction under each condition. Blank cells are missing evidence.

## 4. Input Data Requirements

See Stage 1. Values at reference conditions, uncertainties, and concentrations must come from experiments, curated databases, or a thermodynamic estimator, and they must be consistent in units (kJ/mol, mol/L, K). The script does not read or estimate any value from the internet.

## 5. Expected Computational Output

| Output | Description |
|---|---|
| `reaction_feasibility_by_condition.csv` | $\Delta G$, uncertainty, probability, category, and notes per reaction and condition |
| `pathway_feasibility_by_condition.csv` | Bottleneck, total $\Delta G$, joint feasibility probability, or undetermined status |
| `reaction_sensitivity.csv` | Range of $\Delta G$, temperature and ionic-strength sensitivities, sensitivity flag |
| Condition status table | Which conditions ran or failed in software |
| Pathways robust in all conditions | Count with complete evidence |
| Ranking agreement | Kendall $\tau$ against the baseline |
| `reaction_dG_by_condition.png` | Heat map |

## 6. Scientific Assumptions and Limitations

1. Only thermodynamic feasibility is evaluated. It is distinct from enzyme kinetics, regulation, and actual flux, and a thermodynamically feasible pathway may carry no flux.
2. Metabolite concentrations are held fixed across conditions. They would change in a real cell.
3. Charges and the extended Debye–Hückel approximation are simplifications, and specific ion effects such as Mg$^{2+}$ binding are not modelled.
4. pH and magnesium corrections are not computed. Reaction data are used only at their reference pH.
5. Constant $\Delta H'^{\circ}$ is assumed, so heat-capacity effects are ignored, and it is applied after removing the reference ionic strength although it also depends on ionic strength.
6. Uncertainty covers only the standard deviation of $\Delta G'^{\circ}$, with independent errors.
7. The optional kinetics table (`KINETICS`) is interpolated within its measured temperature range and reported next to the thermodynamic results, but it is not used in any feasibility decision. Enzyme-kinetic constraints, FBA, and global sensitivity analysis (for example with SALib) are not implemented.
8. Pathways are enumerated from stoichiometric connectivity only and may include biochemically unrealistic routes. The cap `MAX_PATHS` may hide some routes.

## 7. Conclusion

The script recalculates reaction thermodynamics under user-defined temperature and ionic conditions, using ionic strength rather than salt concentration, and reports which reactions and pathways stay feasible and keep their ranking. Results are limited by the supplied data and describe thermodynamic feasibility only.
