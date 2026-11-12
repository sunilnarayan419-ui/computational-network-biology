from collections import Counter
from itertools import combinations

import cobra
import numpy as np
import pandas as pd
from scipy import optimize
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Configuration (all biochemical data are user supplied; nothing is built in)
STOICH = "stoichiometry.csv"      # reaction_id, metabolite, coefficient (negative = substrate)
REACTIONS = "reactions.csv"       # reaction_id, vmax_M_per_h [, reversible]
RESPONSE = "response.csv"         # reaction_id, variable (pH | T_K | I_M), value, relative_activity (measured)
MEDIUM = "medium.csv"             # species, kind (acid|neutral|strong_cation|strong_anion), conc_M, charge, pka
SCENARIOS = "scenarios.csv"       # scenario_id, kind, magnitude [, mag_low, mag_high, t_start_h, t_end_h, target, weight]
MITIGATIONS = "mitigations.csv"   # mitigation_id, kind, value [, target], cost
SOURCE, TARGET = "pyruvate", "lactate"
EXCHANGE = []                     # other metabolites allowed to cross the system boundary freely (e.g. CO2)
BASE_T_K = 298.15
PKA_T_K = 298.15                  # temperature at which the supplied pKa values and PKW are valid
PKW = 14.0                        # at PKA_T_K
HORIZON_H, N_STEPS = 24.0, 48
MAX_DAMAGE = 0.20                 # acceptable worst-case performance loss
COST_WEIGHT = 1.0                 # objective = damage + COST_WEIGHT * cost
MAX_COMBO = 2                     # largest number of mitigations combined
VULN_FACTOR = 0.5                 # enzyme-activity scale used for the per-reaction vulnerability screen
REL_MAX = 1.5                     # response values above this are treated as invalid by data validation
RECOVERY_TOL = 0.05
MEASURED_INITIAL_PH = None        # optional check of the medium specification
B_DH = 1.6

DT = HORIZON_H / N_STEPS
out_of_range = Counter()

# 1. Load and validate user data
st = pd.read_csv(STOICH).groupby(["reaction_id", "metabolite"], as_index=False)["coefficient"].sum()
rxn = pd.read_csv(REACTIONS).set_index("reaction_id")
resp = pd.read_csv(RESPONSE)
med = pd.read_csv(MEDIUM)
scn = pd.read_csv(SCENARIOS)
mit = pd.read_csv(MITIGATIONS)

if not set(st["reaction_id"]) <= set(rxn.index) or "vmax_M_per_h" not in rxn:
    raise ValueError("Every reaction needs a vmax_M_per_h entry.")
if med.empty or not med["kind"].isin(["acid", "neutral", "strong_cation", "strong_anion"]).all():
    raise ValueError("Medium composition is missing or has an unknown species kind; pH is underdetermined.")
if not {SOURCE, TARGET} <= set(med["species"]) or not {SOURCE, TARGET} <= set(st["metabolite"]):
    raise ValueError("SOURCE and TARGET must appear in both the medium and the stoichiometry.")
if (resp["relative_activity"] < 0).any():
    raise ValueError("Negative relative activity in response data.")

rxn["reversible"] = rxn["reversible"].fillna(False).astype(bool) if "reversible" in rxn else False
for col, default in [("mag_low", np.nan), ("mag_high", np.nan), ("t_start_h", 0.0), ("t_end_h", HORIZON_H),
                     ("target", ""), ("weight", 1.0)]:
    if col not in scn:
        scn[col] = default
for col in ["target", "cost"]:
    if col not in mit:
        mit[col] = "" if col == "target" else 0.0
scn["target"], mit["target"] = scn["target"].fillna(""), mit["target"].fillna("")
scn["mag_low"], scn["mag_high"] = scn["mag_low"].fillna(scn["magnitude"]), scn["mag_high"].fillna(scn["magnitude"])
scn["weight"] = scn["weight"].fillna(1.0) / scn["weight"].fillna(1.0).sum()
if not scn["kind"].isin(["acid_base_pulse", "ionic_salt", "T_shift", "substrate_scale", "enzyme_scale", "data_corrupt"]).all():
    raise ValueError("Unknown scenario kind.")
if not mit["kind"].isin(["ph_control", "T_control", "enzyme_boost", "add_buffer", "data_validation"]).all():
    raise ValueError("Unknown mitigation kind.")

scenarios, mitigations = scn.to_dict("records"), mit.to_dict("records")

rids = list(rxn.index)
curves = {k: (g["value"].to_numpy(float)[np.argsort(g["value"])], g["relative_activity"].to_numpy(float)[np.argsort(g["value"])])
          for k, g in resp.groupby(["reaction_id", "variable"])}

acids = {r.species: (r.charge, [float(x) for x in str(r.pka).split(";")]) for r in med[med["kind"] == "acid"].itertuples()}
cat = med[med["kind"] == "strong_cation"]
ani = med[med["kind"] == "strong_anion"]
strong_eq = (cat["conc_M"] * cat["charge"]).sum() - (ani["conc_M"] * ani["charge"]).sum()
strong_z2 = (cat["conc_M"] * cat["charge"] ** 2).sum() + (ani["conc_M"] * ani["charge"] ** 2).sum()
C0 = med[med["kind"].isin(["acid", "neutral"])].set_index("species")["conc_M"].to_dict()

# 2. Stoichiometric model (COBRApy) with a limited supply of SOURCE and removal of TARGET
mets = {m: cobra.Metabolite(m) for m in st["metabolite"].unique()}
model = cobra.Model("pipeline")
for rid, g in st.groupby("reaction_id"):
    r = cobra.Reaction(rid)
    r.add_metabolites({mets[m]: c for m, c in zip(g["metabolite"], g["coefficient"])})
    model.add_reactions([r])
for m in set(EXCHANGE) | {SOURCE, TARGET}:
    model.add_boundary(mets[m], type="exchange", reaction_id=f"EX_{m}")
model.objective = f"EX_{TARGET}"
r_obj = [model.reactions.get_by_id(r) for r in rids]
ex_src, ex_tgt = model.reactions.get_by_id(f"EX_{SOURCE}"), model.reactions.get_by_id(f"EX_{TARGET}")
ex_other = [model.reactions.get_by_id(f"EX_{m}") for m in EXCHANGE]
S_med = (st.pivot_table(index="metabolite", columns="reaction_id", values="coefficient", aggfunc="sum")
         .reindex(index=list(C0), columns=rids).fillna(0).to_numpy())


def charge_state(name, pH, A):
    """Fractions and charges of the protonation states of an acid (extended Debye-Hueckel ratio of activity coefficients)."""
    z0, pkas = acids[name]
    j = np.arange(1, len(pkas) + 1)
    ln_ratio = np.log(10) * (pH - np.array(pkas)) - A * (2 * (z0 - j) + 1)
    w = np.r_[0.0, np.cumsum(ln_ratio)]
    w = np.exp(w - w.max())
    return w / w.sum(), z0 - np.arange(len(pkas) + 1)


def speciation(C, T, shift=0.0, salt=0.0):
    """pH and ionic strength from electroneutrality. shift = net strong base added (mol/L, negative = acid); salt = inert 1:1 salt."""
    I_fixed = 0.5 * strong_z2 + salt + 0.5 * abs(shift)
    I = I_fixed

    for _ in range(100):
        s = np.sqrt(I)
        A = (1.10708 - 1.54508e-3 * T + 5.95584e-6 * T ** 2) * s / (1 + B_DH * s)

        def residual(pH):
            h, oh = 10 ** -pH * np.exp(A), 10 ** (pH - PKW) * np.exp(A)
            return strong_eq + shift + h - oh + sum(C[n] * np.dot(*charge_state(n, pH, A)) for n in acids)

        if residual(0.0) * residual(14.0) > 0:
            raise ValueError("No electroneutral solution between pH 0 and 14; medium composition insufficient.")

        pH = optimize.brentq(residual, 0.0, 14.0)
        h, oh = 10 ** -pH * np.exp(A), 10 ** (pH - PKW) * np.exp(A)
        I_new = I_fixed + 0.5 * (h + oh) + 0.5 * sum(
            C[n] * np.dot(f, z ** 2) for n in acids for f, z in [charge_state(n, pH, A)])
        if abs(I_new - I) < 1e-10:
            break
        I = I_new

    return pH, I


def simulate(sc=None, mag=None, mits=()):
    """Time-stepped batch run: speciation -> condition-dependent capacities -> FBA -> medium update."""
    C, cur, detected = dict(C0), curves, False
    kind, t0, t1, tgt = (sc["kind"], sc["t_start_h"], sc["t_end_h"], sc["target"]) if sc else (None, 0, 0, "")

    for m in mits:
        if m["kind"] == "add_buffer":
            C[m["target"]] += m["value"]

    if kind == "data_corrupt":
        cur = {k: (xs, ys * mag) if (not tgt or k[0] == tgt) else (xs, ys) for k, (xs, ys) in curves.items()}
        detected = any((ys < 0).any() or (ys > REL_MAX).any() for k, (xs, ys) in cur.items() if not tgt or k[0] == tgt)
        if detected and any(m["kind"] == "data_validation" for m in mits):
            cur = curves                                   # invalid input rejected, original data used

    t_grid, rate, pHs = np.arange(N_STEPS) * DT, np.zeros(N_STEPS), np.zeros(N_STEPS)
    for i, t in enumerate(t_grid):
        active = t0 <= t < t1
        T, shift, salt, src_scale, enz = BASE_T_K, 0.0, 0.0, 1.0, {}
        if kind == "T_shift" and active:
            T += mag
        elif kind == "acid_base_pulse" and t >= t0:
            shift = mag
        elif kind == "ionic_salt" and t >= t0:
            salt = mag
        elif kind == "substrate_scale" and active:
            src_scale = mag
        elif kind == "enzyme_scale" and active:
            enz[tgt or "*"] = mag

        pH_set = None
        for m in mits:
            if m["kind"] == "T_control":
                T = m["value"]
            elif m["kind"] == "ph_control":
                pH_set = m["value"]
            elif m["kind"] == "enzyme_boost":
                enz[m["target"] or "*"] = enz.get(m["target"] or "*", 1.0) * m["value"]

        pH, I = speciation(C, T, shift, salt)
        pH = pH if pH_set is None else pH_set
        pHs[i] = pH

        for r, rid in zip(r_obj, rids):
            f = 1.0
            for var, val in (("pH", pH), ("T_K", T), ("I_M", I)):
                if (rid, var) in cur:
                    xs, ys = cur[(rid, var)]
                    out_of_range[(rid, var)] += not xs[0] <= val <= xs[-1]
                    f *= np.interp(val, xs, ys)
            cap = rxn.loc[rid, "vmax_M_per_h"] * max(f, 0.0) * enz.get(rid, 1.0) * enz.get("*", 1.0)
            r.bounds = (-cap if rxn.loc[rid, "reversible"] else 0.0, cap)

        ex_src.bounds = (-C[SOURCE] * src_scale / DT, 0.0)
        ex_tgt.bounds = (0.0, 1000.0)
        for r in ex_other:
            r.bounds = (-1000.0, 1000.0)

        sol = model.optimize()
        if sol.status != "optimal":
            continue
        rate[i] = sol.fluxes[ex_tgt.id]
        v = sol.fluxes.reindex(rids).fillna(0.0).to_numpy()
        for name, d in zip(C0, S_med @ v * DT):
            C[name] = max(C[name] + d, 0.0)

    return {"t": t_grid, "rate": rate, "pH": pHs, "J": C[TARGET] - C0[TARGET], "detected": detected}


# 3. Baseline medium and performance
pH0, I0 = speciation(C0, BASE_T_K)
beta = 2e-5 / (speciation(C0, BASE_T_K, 1e-5)[0] - speciation(C0, BASE_T_K, -1e-5)[0])
base = simulate()
J0 = base["J"]

print(f"Baseline medium: pH = {pH0:.2f}, ionic strength = {I0:.4f} M, buffer capacity = {beta:.2e} mol/L per pH unit")
if MEASURED_INITIAL_PH is not None and abs(MEASURED_INITIAL_PH - pH0) > 0.3:
    print(f"WARNING: modelled initial pH differs from the measured value {MEASURED_INITIAL_PH}; composition may be incomplete.")
if abs(BASE_T_K - PKA_T_K) > 5:
    print("WARNING: pKa values and PKW are used unchanged away from PKA_T_K.")
if J0 <= 1e-12:
    raise ValueError("Baseline performance (product formed) is zero; damage is undefined. Check stoichiometry and capacities.")
print(f"Baseline performance J0 = {J0:.4g} mol/L of {TARGET} formed in {HORIZON_H:g} h")

unassessed = [r for r in rids if not any((r, v) in curves for v in ("pH", "T_K", "I_M"))]
if unassessed:
    print("No response data for", unassessed, "- not assumed insensitive; unassessed for pH, temperature and ionic strength.")

# 4. Per-reaction vulnerability: performance loss when one enzyme's activity is scaled by VULN_FACTOR
vuln = []
for rid in rids:
    sc = {"kind": "enzyme_scale", "t_start_h": 0.0, "t_end_h": HORIZON_H, "target": rid}
    vuln.append({"reaction": rid, "damage": 1 - simulate(sc, VULN_FACTOR)["J"] / J0, "has_response_data": rid not in unassessed})
vuln = pd.DataFrame(vuln).sort_values("damage", ascending=False)
vuln.to_csv("vulnerability_screen.csv", index=False)
print(f"\nReaction vulnerability (damage D = 1 - J/J0 when activity x{VULN_FACTOR}):\n", vuln.round(3).to_string(index=False))

# 5. Scenarios x mitigation strategies, each at low / nominal / high magnitude
combos = [()] + [c for k in range(1, MAX_COMBO + 1) for c in combinations(range(len(mitigations)), k)]
rows, failures, traj = [], [], {}

for sc in scenarios:
    cls = "data_integrity" if sc["kind"] == "data_corrupt" else "biochemical"
    start = sc["t_start_h"] if sc["kind"] in ("acid_base_pulse", "ionic_salt", "data_corrupt") else sc["t_end_h"]
    for combo in combos:
        chosen = [mitigations[i] for i in combo]
        row = {"scenario": sc["scenario_id"], "class": cls, "weight": sc["weight"],
               "mitigation": "+".join(m["mitigation_id"] for m in chosen) or "none", "cost": sum(m["cost"] for m in chosen)}
        try:
            for label, mag in (("low", sc["mag_low"]), ("nominal", sc["magnitude"]), ("high", sc["mag_high"])):
                sim = simulate(sc, mag, chosen)
                row[f"damage_{label}"] = 1 - sim["J"] / J0
                if label == "nominal":
                    traj[(sc["scenario_id"], row["mitigation"])] = sim
                    tail = slice(-max(N_STEPS // 10, 1), None)
                    row["residual_damage"] = 1 - sim["rate"][tail].mean() / max(base["rate"][tail].mean(), 1e-12)
                    ok = (sim["t"] >= start) & (sim["rate"] >= (1 - RECOVERY_TOL) * base["rate"])
                    row["recovery_time_h"] = sim["t"][ok][0] - start if ok.any() else np.nan
                    row["corruption_detected"] = sim["detected"]
        except Exception as err:                          # software/pipeline failure, not biochemical damage
            failures.append((sc["scenario_id"], row["mitigation"], str(err)))
            continue
        row["worst_damage"] = max(row["damage_low"], row["damage_nominal"], row["damage_high"])
        row["objective"] = row["damage_nominal"] + COST_WEIGHT * row["cost"]
        rows.append(row)

res = pd.DataFrame(rows)
if res.empty:
    raise SystemExit(f"All simulations failed: {failures[:3]}")
res.to_csv("scenario_mitigation_results.csv", index=False)
print(f"\nSoftware failures (separate from biochemical damage): {len(failures)}")
for f in failures[:10]:
    print("  ", f)

# 6. Report: unmitigated vs best mitigation per scenario, then strategies across all scenarios
best = res.loc[res.groupby("scenario")["objective"].idxmin()]
none = res[res["mitigation"] == "none"].set_index("scenario")
summary = pd.DataFrame({"class": none["class"], "damage_unmitigated": none["damage_nominal"],
                        "range_unmitigated": none["damage_low"].round(3).astype(str) + " to " + none["damage_high"].round(3).astype(str),
                        "recovery_h_unmitigated": none["recovery_time_h"]})
b = best.set_index("scenario")
summary = summary.join(b[["mitigation", "cost", "damage_nominal", "worst_damage", "residual_damage", "recovery_time_h"]],
                       rsuffix="_best")
print("\nBaseline vs disturbed vs best mitigated (damage D = 1 - J/J0; negative = gain):\n", summary.round(3).to_string())

strategies = pd.DataFrame({name: {
    "cost": g["cost"].iloc[0], "expected_damage": np.average(g["damage_nominal"], weights=g["weight"]),
    "worst_case_damage": g["worst_damage"].max(), "n_scenarios": len(g)} for name, g in res.groupby("mitigation")}).T.astype(float)
strategies["objective"] = strategies["expected_damage"] + COST_WEIGHT * strategies["cost"]
strategies = strategies[strategies["n_scenarios"] == len(scenarios)].sort_values("objective")
strategies.to_csv("strategy_summary.csv")
print("\nStrategies across all scenarios:\n", strategies.round(3).to_string())

feasible = strategies[strategies["worst_case_damage"] <= MAX_DAMAGE]
print(f"\nLowest-objective strategy: {strategies.index[0]}")
print(f"Cheapest strategy keeping worst-case damage <= {MAX_DAMAGE}:",
      feasible.sort_values("cost").index[0] if len(feasible) else "none found")
if out_of_range:
    print("Response curves evaluated outside their measured range (edge value used):", dict(out_of_range))
print("Computational predictions only; they are not experimental proof. Mitigations are evaluated, never prescribed.")

# 7. Figure: damage per scenario and pH trajectories
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
x = np.arange(len(summary))
axes[0].bar(x - 0.2, summary["damage_unmitigated"], 0.4, label="unmitigated")
axes[0].bar(x + 0.2, summary["damage_nominal"], 0.4, label="best mitigation")
axes[0].set_xticks(x, summary.index, rotation=45, ha="right")
axes[0].set(ylabel="performance damage D")
axes[0].legend()

worst = summary["damage_unmitigated"].idxmax()
axes[1].plot(base["t"], base["pH"], label="baseline")
axes[1].plot(traj[(worst, "none")]["t"], traj[(worst, "none")]["pH"], label=f"{worst} (unmitigated)")
axes[1].plot(traj[(worst, b.loc[worst, "mitigation"])]["t"], traj[(worst, b.loc[worst, "mitigation"])]["pH"],
             "--", label=f"{worst} ({b.loc[worst, 'mitigation']})")
axes[1].set(xlabel="time (h)", ylabel="pH")
axes[1].legend(fontsize=7)
plt.tight_layout()
plt.savefig("pipeline_resilience.png", dpi=200)
