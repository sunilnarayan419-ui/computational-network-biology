from itertools import islice

import numpy as np
import pandas as pd
import networkx as nx
from scipy import stats
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Configuration (all biochemical data are user supplied; nothing is built in)
REACTIONS = "reactions.csv"       # reaction_id, dG0_prime_kJ_mol, T_ref_K, I_ref_M, pH_ref [, dG0_sd_kJ_mol, dH0_kJ_mol]
STOICH = "stoichiometry.csv"      # reaction_id, metabolite, coefficient (negative = substrate)
METABOLITES = "metabolites.csv"   # metabolite, charge, conc_M
CONDITIONS = "conditions.csv"     # condition_id, T_K, pH  (first row = baseline)
IONS = "ions.csv"                 # condition_id, ion, conc_M, charge
KINETICS = None                   # optional: reaction_id, T_K, relative_activity (measured)
SOURCE, TARGET = "pyruvate", "lactate"
CURRENCY = []                     # cofactors not used as path links (e.g. ATP, NADH); still used in dG
UNIT_ACTIVITY = []                # species with activity fixed at 1 in the transformed dG0' (e.g. water, H+)
MAX_PATH_LEN, MAX_PATHS = 6, 1000
P_THRESHOLD = 0.95                # probability needed to call a direction feasible
N_MC = 2000                       # Monte Carlo draws for dG0' uncertainty
PH_TOL = 0.05                     # reaction data are used only at their reference pH
SENS_DT, SENS_DI, SENSITIVE_KJ = 10.0, 0.05, 1.0
SEED = 42

R = 8.314462618e-3                # kJ mol^-1 K^-1
B_DH = 1.6                        # L^0.5 mol^-0.5 (Alberty extended Debye-Hueckel)
rng = np.random.default_rng(SEED)

# 1. Load and validate user data
rx = pd.read_csv(REACTIONS).set_index("reaction_id")
st = pd.read_csv(STOICH)
met = pd.read_csv(METABOLITES).set_index("metabolite")
cond = pd.read_csv(CONDITIONS).set_index("condition_id")
ions = pd.read_csv(IONS)

for df, cols in [(rx, ["dG0_prime_kJ_mol", "T_ref_K", "I_ref_M", "pH_ref"]), (st, ["reaction_id", "metabolite", "coefficient"]),
                 (met, ["charge", "conc_M"]), (cond, ["T_K", "pH"]), (ions, ["condition_id", "ion", "conc_M", "charge"])]:
    absent = set(cols) - set(df.columns) - {df.index.name}
    if absent:
        raise ValueError(f"Missing columns {sorted(absent)}")

for col in ["dG0_sd_kJ_mol", "dH0_kJ_mol"]:
    if col not in rx:
        rx[col] = np.nan

st = st.groupby(["reaction_id", "metabolite"], as_index=False)["coefficient"].sum()
st = st[st["coefficient"] != 0]

if not set(st["reaction_id"]) <= set(rx.index):
    raise ValueError("Stoichiometry contains reactions missing from the reaction table.")
if not {SOURCE, TARGET} <= set(st["metabolite"]):
    raise ValueError("SOURCE and TARGET must appear in the stoichiometry.")

nu = {rid: g.set_index("metabolite")["coefficient"]
      for rid, g in st[~st["metabolite"].isin(UNIT_ACTIVITY)].groupby("reaction_id")}
kin = pd.read_csv(KINETICS) if KINETICS else None


def ln_gamma(z, I, T):
    """Extended Debye-Hueckel activity coefficient (Alberty), I in mol/L, T in K."""
    alpha = 1.10708 - 1.54508e-3 * T + 5.95584e-6 * T ** 2
    return -alpha * z ** 2 * np.sqrt(I) / (1 + B_DH * np.sqrt(I))


def delta_g(rid, T, I, pH):
    """Transformed dG of a reaction as written (kJ/mol) at T (K), ionic strength I (M), pH. Returns (dG, flags)."""
    r, v, flags = rx.loc[rid], nu.get(rid, pd.Series(dtype=float)), []

    if pd.isna(r["dG0_prime_kJ_mol"]):
        return np.nan, "missing dG0'"
    if abs(pH - r["pH_ref"]) > PH_TOL:
        return np.nan, "no pH-correction data (pH differs from reference)"

    conc = met["conc_M"].reindex(v.index)
    if conc.isna().any() or (conc <= 0).any():
        return np.nan, "missing concentration"

    z = met["charge"].reindex(v.index)
    if z.isna().any():
        lng_ref = lng = np.zeros(len(v))
        flags.append("ideal solution (charge missing)")
    else:
        lng_ref = ln_gamma(z.to_numpy(), r["I_ref_M"], r["T_ref_K"])
        lng = ln_gamma(z.to_numpy(), I, T)

    Tref = r["T_ref_K"]
    g0 = r["dG0_prime_kJ_mol"] - R * Tref * np.dot(v, lng_ref)       # remove reference ionic-strength effect
    if pd.notna(r["dH0_kJ_mol"]):
        g0 = g0 * T / Tref + r["dH0_kJ_mol"] * (1 - T / Tref)        # Gibbs-Helmholtz, constant dH
    elif abs(T - Tref) > 1e-9:
        flags.append("dH missing: only RT ln Q updated")

    return g0 + R * T * np.dot(v, lng + np.log(conc.to_numpy())), "; ".join(flags)


# 2. Candidate pathways SOURCE -> TARGET (either direction of each reaction; currency/unit species not used as links)
Gp = nx.DiGraph()
for rid, v in nu.items():
    v = v[~v.index.isin(CURRENCY)]
    for d in (1, -1):
        node = (rid, d)
        Gp.add_edges_from((m, node) for m, c in v.items() if c * d < 0)
        Gp.add_edges_from((node, m) for m, c in v.items() if c * d > 0)

paths = []
if SOURCE in Gp and TARGET in Gp:
    for p in islice(nx.all_simple_paths(Gp, SOURCE, TARGET, cutoff=2 * MAX_PATH_LEN), MAX_PATHS):
        steps = [n for n in p if isinstance(n, tuple)]
        if len({s[0] for s in steps}) == len(steps):
            paths.append(steps)

if not paths:
    raise ValueError(f"No pathway from {SOURCE} to {TARGET} within {MAX_PATH_LEN} reactions.")

path_rxns = sorted({rid for p in paths for rid, _ in p})
col = {rid: i for i, rid in enumerate(path_rxns)}
eps = rng.standard_normal((N_MC, len(path_rxns)))          # common random numbers across conditions
print(f"{len(paths)} candidate pathways over {len(path_rxns)} reactions (capped at {MAX_PATHS}).")

# 3. Ionic strength and thermodynamics under each condition (failures are isolated per condition)
rxn_rows, path_rows, status = [], [], []

for cid, c in cond.iterrows():
    try:
        T, pH = float(c["T_K"]), float(c["pH"])
        if not (T > 0 and np.isfinite(pH)):
            raise ValueError("temperature must be > 0 K and pH finite")
        ion = ions[ions["condition_id"] == cid]
        I = 0.5 * (ion["conc_M"] * ion["charge"] ** 2).sum()    # ionic strength from composition, not salt concentration

        dG, sd, flag = {}, {}, {}
        for rid in path_rxns:
            dG[rid], flag[rid] = delta_g(rid, T, I, pH)
            sd[rid] = rx.loc[rid, "dG0_sd_kJ_mol"]
            p_fwd = stats.norm.cdf(-dG[rid] / sd[rid]) if pd.notna(sd[rid]) and sd[rid] > 0 else (np.nan if np.isnan(dG[rid]) else float(dG[rid] < 0))
            category = ("missing evidence" if np.isnan(dG[rid]) else "forward" if p_fwd >= P_THRESHOLD
                        else "reverse" if p_fwd <= 1 - P_THRESHOLD else "near equilibrium / uncertain")
            act = np.nan
            if kin is not None:
                k = kin[kin["reaction_id"] == rid].sort_values("T_K")
                if len(k) > 1 and k["T_K"].min() <= T <= k["T_K"].max():
                    act = np.interp(T, k["T_K"], k["relative_activity"])
            rxn_rows.append({"condition": cid, "T_K": T, "pH": pH, "I_M": I, "reaction": rid, "dG_kJ_mol": dG[rid],
                             "dG0_sd_kJ_mol": sd[rid], "P_forward": p_fwd, "category": category,
                             "relative_activity_measured": act, "notes": flag[rid]})

        for n, steps in enumerate(paths):
            d_steps = np.array([d * dG[rid] for rid, d in steps])
            if np.isnan(d_steps).any():
                path_rows.append({"condition": cid, "pathway": n, "status": "undetermined (missing evidence)"})
                continue
            sd_steps = np.nan_to_num([rx.loc[rid, "dG0_sd_kJ_mol"] for rid, _ in steps])
            draws = d_steps + np.array([d for _, d in steps]) * sd_steps * eps[:, [col[rid] for rid, _ in steps]]
            path_rows.append({"condition": cid, "pathway": n, "status": "evaluated",
                              "reactions": " > ".join(f"{rid}{'' if d > 0 else '(rev)'}" for rid, d in steps),
                              "bottleneck_dG_kJ_mol": d_steps.max(), "total_dG_kJ_mol": d_steps.sum(),
                              "P_all_steps_feasible": (draws < 0).all(axis=1).mean()})
        status.append({"condition": cid, "status": "ok", "I_M": I, "message": ""})
    except Exception as err:
        status.append({"condition": cid, "status": "failed", "I_M": np.nan, "message": str(err)})

rxn_df, path_df, status_df = pd.DataFrame(rxn_rows), pd.DataFrame(path_rows), pd.DataFrame(status)
rxn_df.to_csv("reaction_feasibility_by_condition.csv", index=False)
path_df.to_csv("pathway_feasibility_by_condition.csv", index=False)
print("\nCondition status (software failures are separate from biochemical infeasibility):")
print(status_df.to_string(index=False))

if rxn_df.empty:
    raise SystemExit("All conditions failed; see messages above.")

# 4. Reaction robustness and one-factor-at-a-time sensitivity at the baseline condition
base = cond.iloc[0]
I_base = 0.5 * (ions[ions["condition_id"] == cond.index[0]].eval("conc_M * charge ** 2")).sum()
sens = []
for rid in path_rxns:
    g = rxn_df[rxn_df["reaction"] == rid]
    dT = delta_g(rid, base["T_K"] + SENS_DT, I_base, base["pH"])[0] - delta_g(rid, base["T_K"] - SENS_DT, I_base, base["pH"])[0]
    dI = delta_g(rid, base["T_K"], I_base + SENS_DI, base["pH"])[0] - delta_g(rid, base["T_K"], max(I_base - SENS_DI, 0), base["pH"])[0]
    sens.append({"reaction": rid, "dG_range_kJ_mol": g["dG_kJ_mol"].max() - g["dG_kJ_mol"].min(),
                 "categories_seen": ", ".join(sorted(g["category"].unique())),
                 "ddG_temperature_kJ_mol": dT, "ddG_ionic_strength_kJ_mol": dI,
                 "condition_sensitive": bool(max(abs(dT), abs(dI)) >= SENSITIVE_KJ
                                             or ((g["category"] == "forward").any() and (g["category"] == "reverse").any()))})

sens_df = pd.DataFrame(sens)
sens_df.to_csv("reaction_sensitivity.csv", index=False)
print(f"\nReactions flagged condition-sensitive: {int(sens_df['condition_sensitive'].sum())} of {len(sens_df)}")
print(sens_df.round(3).to_string(index=False))

# 5. Pathway ranking and robustness across conditions
ok = path_df[path_df["status"] == "evaluated"]
if not ok.empty:
    wide = ok.pivot(index="pathway", columns="condition", values="bottleneck_dG_kJ_mol")
    prob = ok.pivot(index="pathway", columns="condition", values="P_all_steps_feasible")
    n_cond = len(cond)
    robust = prob.index[(prob.notna().all(axis=1)) & (prob >= P_THRESHOLD).all(axis=1)] if prob.shape[1] == n_cond else []
    print(f"\nPathways feasible (P >= {P_THRESHOLD}) with complete evidence in all {n_cond} conditions: {len(robust)}")

    first = cond.index[0]
    if first in wide:
        for cid in wide.columns[wide.columns != first]:
            both = wide[[first, cid]].dropna()
            if len(both) > 2:
                tau, p = stats.kendalltau(both[first], both[cid])
                print(f"Ranking agreement (Kendall tau) baseline vs {cid}: {tau:.3f} (p = {p:.3g}, {len(both)} pathways)")

    best = ok.sort_values("bottleneck_dG_kJ_mol").groupby("condition").head(3)
    print("\nTop pathways per condition by bottleneck driving force (more negative is better):")
    print(best[["condition", "pathway", "reactions", "bottleneck_dG_kJ_mol", "P_all_steps_feasible"]].round(3).to_string(index=False))

# 6. Figure: dG of pathway reactions across conditions
heat = rxn_df.pivot(index="reaction", columns="condition", values="dG_kJ_mol")
fig, ax = plt.subplots(figsize=(1.2 * heat.shape[1] + 4, 0.35 * heat.shape[0] + 2))
im = ax.imshow(heat.to_numpy(float), aspect="auto", cmap="coolwarm")
ax.set_xticks(range(heat.shape[1]), heat.columns, rotation=45, ha="right")
ax.set_yticks(range(heat.shape[0]), heat.index)
fig.colorbar(im, label="dG as written (kJ/mol); blank = missing evidence")
plt.tight_layout()
plt.savefig("reaction_dG_by_condition.png", dpi=200)

print("\nThermodynamic feasibility only: enzyme kinetics, regulation and actual flux are not evaluated.")
print("Missing data are reported as missing evidence, never as infeasibility.")
