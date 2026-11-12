import itertools
import re
import time
from collections import Counter, defaultdict

import cobra
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import powerlaw
from scipy.optimize import linprog

# Configuration
MODEL_FILE = "model.xml"                # SBML reconstruction (stoichiometry, bounds, subsystems)
FEAS_CSV = "reaction_feasibility.csv"   # Problem 05 output (reaction_id, direction, admissible, cost); None = bounds only
HIERARCHY_CSV = "hierarchy.csv"         # reaction_id, subsystem, pathway (one row per membership); None = model subsystems
UPDATES_CSV = None                      # action (add|remove), reaction_id [, stoichiometry, lower_bound, upper_bound, subsystem, pathway]
SOURCE, TARGET = "pyr_c", "accoa_c"
POOL = {"h", "h2o", "atp", "adp", "pi", "ppi", "amp", "nad", "nadh", "nadp", "nadph", "coa", "co2"}
UNKNOWN_COST = 0.7        # cost of a bound-allowed direction without feasibility evidence (not "infeasible")
M_OTHER = 0.1             # module-compatibility factor M(i) for metabolites outside the target pathway
K_COARSE = 3              # routes kept per start/end pair at the subsystem and pathway levels
K_PATHS, MAX_LEN = 50, 8
TIE_BREAK = 1e-3
FLUX_CAP = 1000.0
MIN_FIT_NODES = 50        # no power-law fit below this number of nodes
SEED = 42
OUT_PATHS, OUT_FIG = "hierarchical_pathways.csv", "hierarchy_levels.png"


def base(x):
    return re.sub(r"_[a-z]$", "", x.lower())


def allowed(v, d):
    return v["ub"] > 0 if d == "fwd" else v["lb"] < 0


def nonpool(r):
    return {x for x in store[r]["stoich"] if base(x) not in POOL}


# 1. Load reactions (stoichiometry retained) and hierarchy memberships
model = cobra.io.read_sbml_model(MODEL_FILE)
store = {r.id: {"stoich": {x.id: c for x, c in r.metabolites.items()}, "lb": r.lower_bound, "ub": r.upper_bound}
         for r in model.reactions if not r.boundary}
if HIERARCHY_CSV:
    mem = pd.read_csv(HIERARCHY_CSV, dtype=str).fillna("")
    if not {"reaction_id", "subsystem", "pathway"} <= set(mem.columns):
        raise ValueError("Hierarchy file needs: reaction_id, subsystem, pathway.")
else:
    print("No hierarchy file: model subsystems are used and pathway = subsystem (no pathways are invented).")
    mem = pd.DataFrame([(r.id, r.subsystem, r.subsystem) for r in model.reactions if not r.boundary and r.subsystem],
                       columns=["reaction_id", "subsystem", "pathway"])
mem = mem[["reaction_id", "subsystem", "pathway"]]
mem = mem[(mem != "").all(axis=1)]


def check():
    return {"unknown_reactions": len(set(mem["reaction_id"]) - set(store)),
            "duplicate_rows": int(mem.duplicated().sum()),
            "unannotated_reactions": len(set(store) - set(mem["reaction_id"]))}


print("Hierarchy consistency (initial):", check())
mem = mem.drop_duplicates()
mem = mem[mem["reaction_id"].isin(store)]

# 2. Reaction costs C(r) = -ln F(r) from Problem 05; missing evidence gets UNKNOWN_COST, infeasible is excluded
cost = {(r, d): UNKNOWN_COST for r, v in store.items() for d in ("fwd", "rev") if allowed(v, d)}
if FEAS_CSV:
    f = pd.read_csv(FEAS_CSV)
    f["admissible"] = f["admissible"].astype(str).str.lower().eq("true")
    for r in f.itertuples():
        k = (r.reaction_id, r.direction)
        if k in cost:
            if r.admissible:
                cost[k] = float(r.cost)
            else:
                cost.pop(k)
else:
    print("No feasibility table: flux-bound directionality only; all reactions get the same cost.")


def degrees():
    return Counter(x for r in store for x in nonpool(r))


# 3. Dynamic hierarchical updates (hierarchy is derived from reactions + memberships, so it stays consistent)
def snapshot(label):
    return {"update": label, "reactions": len(store), "subsystems": mem["subsystem"].nunique(),
            "pathways": len(mem[["subsystem", "pathway"]].drop_duplicates()), **check()}


log = [snapshot("initial")]
if UPDATES_CSV:
    ups = pd.read_csv(UPDATES_CSV, dtype=str).reindex(columns=[
        "action", "reaction_id", "stoichiometry", "lower_bound", "upper_bound", "subsystem", "pathway"]).fillna("")
    for u in ups.itertuples():
        if u.action == "add":
            if u.reaction_id in store or not (u.stoichiometry and u.subsystem and u.pathway):
                raise ValueError(f"Cannot add reaction {u.reaction_id} (needs stoichiometry, subsystem, pathway).")
            st = {k: float(v) for k, v in (t.rsplit(":", 1) for t in u.stoichiometry.split(";"))}
            # Hierarchical preferential-attachment diagnostic: Pi(i) ~ k(i) B(i) M(i), B = 1 for non-pool metabolites
            deg = degrees()
            members = set(mem.loc[(mem["subsystem"] == u.subsystem) & (mem["pathway"] == u.pathway), "reaction_id"])
            module = {x for r in members for x in nonpool(r)}
            w = {x: k * (1.0 if x in module else M_OTHER) for x, k in deg.items()}
            total = sum(w.values())
            print(f"\nAdd {u.reaction_id} to {u.subsystem} | {u.pathway}: Pi(i) before the update")
            for x in st:
                if base(x) not in POOL:
                    print(f"  {x}: k = {deg.get(x, 0)}, in module = {x in module}, Pi = {w.get(x, 0) / total:.4f}")
            store[u.reaction_id] = {"stoich": st, "lb": float(u.lower_bound), "ub": float(u.upper_bound)}
            mem = pd.concat([mem, pd.DataFrame([{"reaction_id": u.reaction_id, "subsystem": u.subsystem,
                                                 "pathway": u.pathway}])], ignore_index=True)
            for d in ("fwd", "rev"):
                if allowed(store[u.reaction_id], d):
                    cost[(u.reaction_id, d)] = UNKNOWN_COST
        elif u.action == "remove" and u.reaction_id in store:
            store.pop(u.reaction_id)
            mem = mem[mem["reaction_id"] != u.reaction_id]
            for d in ("fwd", "rev"):
                cost.pop((u.reaction_id, d), None)
        else:
            raise ValueError(f"Invalid update: {u.action} {u.reaction_id}")
        log.append(snapshot(f"{u.action} {u.reaction_id}"))
print("\nUpdate log:")
print(pd.DataFrame(log).to_string(index=False))

# 4. Four hierarchical levels with explicit cross-level relationships
mem = mem.assign(pid=mem["subsystem"] + " | " + mem["pathway"])
G4 = nx.Graph()   # level 4/3 backbone: reaction-metabolite bipartite graph (non-pool metabolites)
for r, v in store.items():
    for x in nonpool(r):
        G4.add_edge(("R", r), ("M", x), coef=v["stoich"][x])
G3 = nx.bipartite.projected_graph(G4, [n for n in G4 if n[0] == "R"])   # reactions sharing a metabolite


def group_graph(groups):
    """Groups (pathways/subsystems) are linked only when their reactions share a non-pool metabolite."""
    inv = defaultdict(set)
    for g, rs in groups.items():
        for x in {x for r in rs for x in nonpool(r)}:
            inv[x].add(g)
    G = nx.Graph()
    G.add_nodes_from(groups)
    for x, gs in inv.items():
        for a, b in itertools.combinations(sorted(gs), 2):
            if G.has_edge(a, b):
                G[a][b]["shared"].add(x)
            else:
                G.add_edge(a, b, shared={x})
    return G


pw = {p: set(g["reaction_id"]) for p, g in mem.groupby("pid")}        # overlapping memberships preserved
ss = {s: set(g["reaction_id"]) for s, g in mem.groupby("subsystem")}
sub_of = dict(zip(mem["pid"], mem["subsystem"]))
G1, G2 = group_graph(ss), group_graph(pw)


def describe(name, degs):
    degs = np.asarray(degs)
    pos = degs[degs > 0]
    if len(pos) == 0:
        return print(f"{name}: no connected nodes")
    line = f"{name}: n = {len(degs)}, mean k = {pos.mean():.2f}, max k = {pos.max()}, <k^2>/<k> = {(pos ** 2).mean() / pos.mean():.2f}"
    if len(degs) >= MIN_FIT_NODES and len(set(pos)) >= 5:
        try:
            fit = powerlaw.Fit(pos, discrete=True, verbose=False)
            R, p = fit.distribution_compare("power_law", "lognormal", normalized_ratio=True)
            line += f"; power-law alpha = {fit.power_law.alpha:.2f}, xmin = {fit.xmin}; vs lognormal R = {R:.2f}, p = {p:.3f}"
        except Exception as e:
            line += f"; power-law fit failed ({e})"
    else:
        line += "; too few nodes/distinct degrees for a power-law fit (no scale-free claim)"
    print(line)


print("\nDegree statistics per level (each level assessed independently):")
describe("Level 1 subsystems", [d for _, d in G1.degree()])
describe("Level 2 pathways", [d for _, d in G2.degree()])
describe("Level 3 reactions", [d for _, d in G3.degree()])
describe("Level 4 metabolites (reactions per metabolite)", [d for n, d in G4.degree() if n[0] == "M"])


# 5. Pathway search and stoichiometric check (same definitions as Problem 05)
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


def stoich_check(path):
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


def search(label, reactions=None):
    t0 = time.perf_counter()
    G = search_graph({k: c for k, c in cost.items() if reactions is None or k[0] in reactions})
    s, t = ("M", SOURCE), ("M", TARGET)
    paths = []
    if s in G and t in G:
        try:
            gen = nx.shortest_simple_paths(G, s, t, weight="weight")
            paths = [[(n[1], n[2]) for n in p if n[0] == "R"] for p in itertools.islice(gen, K_PATHS)]
        except nx.NetworkXNoPath:
            pass
    rows = []
    for p in (p for p in paths if len(p) <= MAX_LEN):
        ok, unmet = stoich_check(p)
        mm = mem[mem["reaction_id"].isin([r for r, _ in p])]
        rows.append({"search": label, "pathway": " > ".join(f"{r}({d})" for r, d in p), "n_reactions": len(p),
                     "cost": sum(cost[k] for k in p), "stoichiometrically_supported": ok, "unmet_cosubstrates": unmet,
                     "bottleneck_reaction": max(p, key=lambda k: cost[k])[0],
                     "subsystems": "; ".join(sorted(mm["subsystem"].unique())),
                     "pathways": "; ".join(sorted(mm["pid"].unique()))})
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values(["stoichiometrically_supported", "cost", "n_reactions"], ascending=[False, True, True])
    return df, G.number_of_edges(), time.perf_counter() - t0


def routes(G, starts, ends):
    keep = set()
    for a in starts:
        for b in ends:
            if a == b:
                keep.add(a)
                continue
            try:
                for p in itertools.islice(nx.shortest_simple_paths(G, a, b), K_COARSE):
                    keep.update(p)
            except (nx.NetworkXNoPath, nx.NodeNotFound):
                pass
    return keep


# Coarse-to-fine exploration: subsystems -> pathways -> reactions
src = {r for r in store if SOURCE in store[r]["stoich"]}
tgt = {r for r in store if TARGET in store[r]["stoich"]}
at = lambda rs, col: set(mem.loc[mem["reaction_id"].isin(rs), col])
S_sel = routes(G1, at(src, "subsystem"), at(tgt, "subsystem"))
G2s = G2.subgraph([p for p in G2 if sub_of[p] in S_sel])
P_sel = routes(G2s, at(src, "pid") & set(G2s), at(tgt, "pid") & set(G2s))
allowed_rxns = set().union(*(pw[p] for p in P_sel)) if P_sel else set()
print(f"\nCoarse search: {len(S_sel)} subsystems -> {len(P_sel)} pathways -> {len(allowed_rxns)} of {len(store)} reactions")

hier, e_h, t_h = search("hierarchical", allowed_rxns) if allowed_rxns else (pd.DataFrame(), 0, 0.0)
flat, e_f, t_f = search("non-hierarchical")
out = pd.concat([hier, flat], ignore_index=True)
if len(out):
    out.to_csv(OUT_PATHS, index=False)
best = lambda d: d[d["stoichiometrically_supported"]]["cost"].min() if len(d) else np.nan
print("\nHierarchical vs non-hierarchical search:")
print(pd.DataFrame([
    {"search": "hierarchical", "search_graph_edges": e_h, "seconds": round(t_h, 3), "candidates": len(hier), "best_supported_cost": best(hier)},
    {"search": "non-hierarchical", "search_graph_edges": e_f, "seconds": round(t_f, 3), "candidates": len(flat), "best_supported_cost": best(flat)}]
).to_string(index=False))
if not np.isnan(best(hier)) and not np.isnan(best(flat)):
    print(f"Cost gap (hierarchical - non-hierarchical): {best(hier) - best(flat):.4f}")
for name, d in (("hierarchical", hier), ("non-hierarchical", flat)):
    if len(d):
        print(f"\nTop pathways, {name}:")
        print(d.head(5)[["pathway", "cost", "stoichiometrically_supported", "unmet_cosubstrates", "bottleneck_reaction"]]
              .round(3).to_string(index=False))
hubs = Counter({x: d for (t, x), d in G4.degree() if t == "M"}).most_common(5)
print("\nMost connected metabolites (reactions per metabolite, pool excluded):", hubs)

# 6. Visualization: subsystem and pathway levels, selected modules highlighted
fig, ax = plt.subplots(1, 2, figsize=(15, 6))
for a, G, sel, title in ((ax[0], G1, S_sel, "Level 1: subsystems"), (ax[1], G2s, P_sel, "Level 2: pathways in selected subsystems")):
    nx.draw_networkx(G, pos=nx.spring_layout(G, seed=SEED), ax=a, font_size=6, width=0.3, node_size=250,
                     node_color=["tab:orange" if n in sel else "lightgray" for n in G])
    a.set_title(title + " (selected = orange)")
    a.axis("off")
fig.tight_layout()
fig.savefig(OUT_FIG, dpi=150)
print(f"\nSaved {OUT_FIG}" + (f" and {OUT_PATHS}" if len(out) else ""))
