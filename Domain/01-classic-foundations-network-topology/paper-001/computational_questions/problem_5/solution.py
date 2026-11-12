import copy
import hashlib
import itertools
import re
from collections import Counter

import cobra
import networkx as nx
import numpy as np
import pandas as pd
from cobra.flux_analysis import pfba
from scipy.optimize import linprog

# Configuration
MODEL_FILE = "model.xml"        # SBML reconstruction (stoichiometry, bounds, gene rules)
THERMO_CSV = "thermo.csv"       # reaction_id, dG0_prime_kJ_mol [, dG0_sd_kJ_mol]; None if unavailable
CONC_CSV = "concentrations.csv" # metabolite, conc_min_M, conc_max_M; None if unavailable
UPDATES_CSV = None              # action (add|remove|set_bounds), reaction_id [, stoichiometry, lower_bound, upper_bound]
SOURCE, TARGET = "pyr_c", "accoa_c"
POOL = {"h", "h2o", "atp", "adp", "pi", "ppi", "amp", "nad", "nadh", "nadp", "nadph", "coa", "co2"}
UNAVAILABLE_METABOLITES = set()  # substrates/cofactors absent under the studied condition
UNAVAILABLE_GENES = set()        # genes/enzymes absent under the studied condition
FBA_MEDIUM = {}                  # e.g. {"EX_pyr_e": 10}; empty = keep the model's medium
TEMPERATURE_K = 298.15
Z_SD = 1.96                      # width of the dG0' uncertainty interval (in standard errors)
CONC_DEFAULT = (1e-6, 1e-2)      # M; assumption used only for metabolites missing in CONC_CSV
W_UNKNOWN = 0.5                  # score weight for missing evidence (ordinal, not a probability)
W_MIN = 0.05
TIE_BREAK = 1e-3                 # tiny per-reaction cost so shorter equal-cost pathways rank first
K_PATHS, MAX_LEN = 50, 8
FLUX_CAP = 1000.0
OUT_FEAS, OUT_PATHS = "reaction_feasibility.csv", "ranked_pathways.csv"
R_GAS = 8.314e-3                 # kJ/(mol K)


def base(x):
    return re.sub(r"_[a-z]$", "", x.lower())


def comp(x):
    return x[-1] if re.search(r"_[a-z]$", x) else ""


def digest(s):
    return hashlib.sha256(repr(sorted((r, sorted(v["stoich"].items()), v["lb"], v["ub"])
                                      for r, v in s.items())).encode()).hexdigest()[:16]


def degrees(s):
    return Counter(x for v in s.values() for x in v["stoich"] if base(x) not in POOL)


# 1. Load the reaction store (stoichiometry is retained for every reaction)
model = cobra.io.read_sbml_model(MODEL_FILE)


def enzyme_ok(r):
    if not r.genes:
        return None   # no annotation: missing evidence
    expr = re.sub(r"[^\s()]+", lambda t: t.group().lower() if t.group().lower() in ("and", "or")
                  else str(t.group() not in UNAVAILABLE_GENES), r.gene_reaction_rule)
    return bool(eval(expr, {"__builtins__": {}}))


store = {r.id: {"stoich": {x.id: c for x, c in r.metabolites.items()}, "lb": r.lower_bound,
                "ub": r.upper_bound, "enzyme": enzyme_ok(r)}
         for r in model.reactions if not r.boundary}
if not any(SOURCE in v["stoich"] for v in store.values()) or not any(TARGET in v["stoich"] for v in store.values()):
    raise ValueError("SOURCE or TARGET metabolite not found in the model.")
thermo = pd.read_csv(THERMO_CSV).set_index("reaction_id") if THERMO_CSV else pd.DataFrame()
conc = pd.read_csv(CONC_CSV).set_index("metabolite") if CONC_CSV else pd.DataFrame()


# 2. Reaction feasibility F(r) = f(T, S, E, C, D) per direction
def dg_interval(rid, st):
    """dG' interval (kJ/mol) for the reaction as written; None if dG0' is missing or it is a transport reaction."""
    if rid not in thermo.index or pd.isna(thermo.at[rid, "dG0_prime_kJ_mol"]):
        return None
    if len({comp(x) for x in st}) > 1:
        return None   # transport needs membrane potential / pH gradient, not supplied
    sd = thermo.at[rid, "dG0_sd_kJ_mol"] if "dG0_sd_kJ_mol" in thermo else np.nan
    sd = 0.0 if pd.isna(sd) else sd
    lo, hi = (thermo.at[rid, "dG0_prime_kJ_mol"] + s * Z_SD * sd for s in (-1, 1))
    for x, v in st.items():
        if base(x) in {"h", "h2o"}:
            continue   # dG0' is a transformed value: H+ and water activities are not in Q
        c_lo, c_hi = (conc.at[x, "conc_min_M"], conc.at[x, "conc_max_M"]) if x in conc.index else CONC_DEFAULT
        lo += R_GAS * TEMPERATURE_K * v * np.log(c_lo if v > 0 else c_hi)
        hi += R_GAS * TEMPERATURE_K * v * np.log(c_hi if v > 0 else c_lo)
    return lo, hi


def assess(rid, rec, d):
    sign = 1 if d == "fwd" else -1
    reasons, missing, F, dg_txt = [], [], 1.0, ""
    if (rec["ub"] <= 0) if d == "fwd" else (rec["lb"] >= 0):
        reasons.append("direction blocked by flux bounds")
    lacking = [x for x, c in rec["stoich"].items() if sign * c < 0 and x in UNAVAILABLE_METABOLITES]
    if lacking:
        reasons.append("unavailable substrate/cofactor: " + ",".join(lacking))
    if rec["enzyme"] is False:
        reasons.append("required enzyme/gene unavailable")
    elif rec["enzyme"] is None:
        missing.append("enzyme annotation")
        F *= W_UNKNOWN
    iv = dg_interval(rid, rec["stoich"])
    if iv is None:
        missing.append("dG0'/thermodynamic data")
        F *= W_UNKNOWN
    else:
        lo, hi = iv if sign == 1 else (-iv[1], -iv[0])
        dg_txt = f"{lo:.1f}..{hi:.1f}"
        if lo > 0:
            reasons.append(f"thermodynamically unfavourable (dG' {dg_txt} kJ/mol)")
        elif hi >= 0:   # interval straddles zero: weight = fraction of the interval below zero
            F *= max(W_MIN, -lo / (hi - lo) if hi > lo else 0.5)
    ok = not reasons
    return {"reaction_id": rid, "direction": d, "admissible": ok, "F": F if ok else 0.0,
            "cost": -np.log(F) if ok else np.inf, "dG_range_kJ_mol": dg_txt,
            "reasons": "; ".join(reasons), "missing_evidence": "; ".join(missing)}


def evaluate(s):
    return pd.DataFrame([assess(r, v, d) for r, v in s.items() for d in ("fwd", "rev")])


# 3. Dynamic updates with immutable version history and preferential-attachment diagnostics
history, log = [copy.deepcopy(store)], [("initial", len(store), digest(store))]
if UPDATES_CSV:
    ups = pd.read_csv(UPDATES_CSV, dtype=str).reindex(
        columns=["action", "reaction_id", "stoichiometry", "lower_bound", "upper_bound"]).fillna("")
    for u in ups.itertuples():
        if u.action == "add":
            if u.reaction_id in store or not u.stoichiometry:
                raise ValueError(f"Cannot add reaction {u.reaction_id}.")
            st = {k: float(v) for k, v in (t.rsplit(":", 1) for t in u.stoichiometry.split(";"))}
            deg = degrees(store)
            total = sum(deg.values())   # B(i) = 1 for non-pool metabolites, 0 for pool metabolites
            print(f"\nAdd {u.reaction_id}: attachment propensity Pi(i) = k(i)B(i)/sum(k(j)B(j)) before the update")
            for x in st:
                if base(x) not in POOL:
                    print(f"  {x}: k = {deg.get(x, 0)}, Pi = {deg.get(x, 0) / total:.4f}")
            store[u.reaction_id] = {"stoich": st, "lb": float(u.lower_bound), "ub": float(u.upper_bound),
                                    "enzyme": None}
        elif u.reaction_id not in store:
            raise ValueError(f"Reaction {u.reaction_id} not in network.")
        elif u.action == "remove":
            store.pop(u.reaction_id)
        elif u.action == "set_bounds":
            store[u.reaction_id]["lb"], store[u.reaction_id]["ub"] = float(u.lower_bound), float(u.upper_bound)
        else:
            raise ValueError(f"Unknown action {u.action}.")
        history.append(copy.deepcopy(store))
        log.append((f"{u.action} {u.reaction_id}", len(store), digest(store)))
history_intact = all(digest(h) == l[2] for h, l in zip(history, log))
print("\nNetwork versions:")
for i, l in enumerate(log):
    print(f"  v{i}: {l[0]} | {l[1]} reactions | digest {l[2]}")
print("Historical versions unchanged:", history_intact)
deg = degrees(store)
print("Top metabolites by reaction participation (pool excluded):", deg.most_common(5))

# 4. Evaluate feasibility of the final network version
feas = evaluate(store)
feas.to_csv(OUT_FEAS, index=False)
rejected = feas[~feas["admissible"] & (feas["reasons"] != "direction blocked by flux bounds")]
print(f"\nReaction directions: {len(feas)} | admissible: {int(feas['admissible'].sum())} | "
      f"rejected with reasons: {len(rejected)} | admissible with missing evidence: "
      f"{int((feas['admissible'] & (feas['missing_evidence'] != '')).sum())}")
print(rejected.head(10)[["reaction_id", "direction", "reasons"]].to_string(index=False))


# 5. Pathway search on the bipartite metabolite-reaction graph (reaction identity is kept)
def search_graph(weights):
    G = nx.DiGraph()
    for (rid, d), c in weights.items():
        sign = 1 if d == "fwd" else -1
        st = store[rid]["stoich"]
        subs = [x for x, v in st.items() if sign * v < 0 and base(x) not in POOL]
        prods = [x for x, v in st.items() if sign * v > 0 and base(x) not in POOL]
        if subs and prods:
            for x in subs:
                G.add_edge(("M", x), ("R", rid, d), weight=c + TIE_BREAK)
            for x in prods:
                G.add_edge(("R", rid, d), ("M", x), weight=0.0)
    return G


def top_paths(G, k):
    s, t = ("M", SOURCE), ("M", TARGET)
    if s not in G or t not in G:
        return []
    try:
        paths = list(itertools.islice(nx.shortest_simple_paths(G, s, t, weight="weight"), k))
    except nx.NetworkXNoPath:
        return []
    return [[(n[1], n[2]) for n in p if n[0] == "R"] for p in paths]


def stoich_check(path):
    """LP: can these oriented reactions run at steady state and make 1 unit of TARGET?
    Intermediates must not be net consumed; pool metabolites and SOURCE are assumed supplied."""
    sign = lambda d: 1 if d == "fwd" else -1
    mets = sorted({x for r, d in path for x in store[r]["stoich"]})
    free = lambda x: base(x) in POOL or x == SOURCE
    S = np.array([[sign(d) * store[r]["stoich"].get(x, 0.0) for r, d in path] for x in mets])
    rows = [i for i, x in enumerate(mets) if not free(x)]
    bounds = [(0, min(FLUX_CAP, store[r]["ub"] if d == "fwd" else -store[r]["lb"])) for r, d in path]
    res = linprog(np.ones(len(path)), A_ub=-S[rows], b_ub=[-1.0 if mets[i] == TARGET else 0.0 for i in rows],
                  bounds=bounds, method="highs")
    if res.success:
        return True, ""
    made = {x for r, d in path for x, c in store[r]["stoich"].items() if sign(d) * c > 0}
    used = {x for r, d in path for x, c in store[r]["stoich"].items() if sign(d) * c < 0}
    return False, ",".join(sorted(x for x in used - made if not free(x))) or "no balanced steady state"


# FBA reference on the unmodified input model (not on dynamic updates)
fba_flux, fba_max = {}, None
with model:
    if FBA_MEDIUM:
        model.medium = FBA_MEDIUM
    model.objective = model.add_boundary(model.metabolites.get_by_id(TARGET), type="demand")
    sol = model.optimize()
    if sol.status == "optimal" and sol.objective_value > 1e-9:
        fba_max, fba_flux = sol.objective_value, pfba(model).fluxes.to_dict()
print("\nFBA maximum TARGET production under the chosen medium:", fba_max if fba_max else "none / infeasible")

idx = feas.set_index(["reaction_id", "direction"])
weights = {k: c for k, c in idx["cost"].items() if np.isfinite(c)}
results = []
for path in (p for p in top_paths(search_graph(weights), K_PATHS) if len(p) <= MAX_LEN):
    ok, unmet = stoich_check(path)
    rows = [idx.loc[k] for k in path]
    active = [fba_flux.get(r, 0.0) * (1 if d == "fwd" else -1) > 1e-9 for r, d in path] if fba_flux else None
    results.append({
        "pathway": " > ".join(f"{r}({d})" for r, d in path), "n_reactions": len(path),
        "cost": sum(x["cost"] for x in rows), "stoichiometrically_supported": ok, "unmet_cosubstrates": unmet,
        "missing_evidence": "; ".join(sorted({e for x in rows for e in x["missing_evidence"].split("; ") if e})),
        "dG_ranges_kJ_mol": " | ".join(x["dG_range_kJ_mol"] or "NA" for x in rows),
        "fba_active_reactions": f"{sum(active)}/{len(path)}" if active else "NA"})
if results:
    out = pd.DataFrame(results).sort_values(["stoichiometrically_supported", "cost", "n_reactions"],
                                            ascending=[False, True, True])
    out.insert(0, "rank", range(1, len(out) + 1))
    out.to_csv(OUT_PATHS, index=False)
    print(f"\nCandidate pathways {SOURCE} -> {TARGET} (supported pathways ranked first):")
    print(out.head(10).round(3).to_string(index=False))
else:
    print(f"\nNo admissible pathway from {SOURCE} to {TARGET} within {MAX_LEN} reactions.")

# 6. Baseline: unweighted shortest path using flux bounds only (no thermodynamic/enzyme/availability checks)
allowed = {(r, d): 1.0 for r, v in store.items() for d in ("fwd", "rev") if (v["ub"] > 0 if d == "fwd" else v["lb"] < 0)}
base_paths = top_paths(search_graph(allowed), 1)
if base_paths:
    bp = base_paths[0]
    print("\nUnconstrained shortest path:", " > ".join(f"{r}({d})" for r, d in bp))
    for k in bp:
        print(f"  {k[0]}({k[1]}): admissible = {idx.loc[k, 'admissible']}; {idx.loc[k, 'reasons'] or 'no rejection reason'}")
    print("  stoichiometrically supported:", stoich_check(bp))
print(f"\nSaved {OUT_FEAS}" + (f" and {OUT_PATHS}" if results else ""))
