import hashlib
import json
import re

import networkx as nx
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.model_selection import GroupShuffleSplit
from sklearn.preprocessing import StandardScaler

# Configuration
INPUT = "expression_long.csv"   # columns: gene, sample, time, expression [, group]
AUDIT_OUTPUT = "leakage_audit.json"
TRAIN_FRAC, VAL_FRAC = 0.6, 0.2  # test fraction = remainder
TOP_GENES = 500                  # variance-based gene selection (training data only)
MIN_ABS_CORR = 0.7               # fixed a priori; never tuned on validation/test data
MIN_SAMPLES = 5                  # minimum samples needed to estimate correlations
SEED = 42

# 1. Load and validate the long-format expression table
df = pd.read_csv(INPUT, dtype={"gene": str, "sample": str, "time": str, "group": str})
if not {"gene", "sample", "time", "expression"} <= set(df.columns):
    raise ValueError("Required columns: gene, sample, time, expression.")
if df[["gene", "sample", "time"]].isna().any().any():
    raise ValueError("gene, sample and time must not be missing.")
if df.duplicated(["gene", "sample"]).any():
    raise ValueError("Duplicate gene-sample records; resolve replicates first.")
if TRAIN_FRAC <= 0 or VAL_FRAC <= 0 or TRAIN_FRAC + VAL_FRAC >= 1:
    raise ValueError("Need TRAIN_FRAC > 0, VAL_FRAC > 0 and TRAIN_FRAC + VAL_FRAC < 1.")

df["expression"] = pd.to_numeric(df["expression"], errors="coerce")
if "group" not in df:
    df["group"] = df["sample"]  # no subject/replicate info: each sample is its own group
if (df.groupby("sample")[["time", "group"]].nunique() > 1).to_numpy().any():
    raise ValueError("Each sample must have exactly one time point and one group.")

meta = df.drop_duplicates("sample").set_index("sample")[["time", "group"]]
X = df.pivot(index="sample", columns="gene", values="expression").loc[meta.index]


def natural(s):
    return [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", s)]


times = sorted(meta["time"].unique(), key=natural)
meta["t"] = meta["time"].map({t: i for i, t in enumerate(times)})
t_of = meta["t"]

# 2. Leakage-controlled partitioning
rank = {"train": 0, "val": 1, "test": 2}
temporal = len(times) >= 3
if temporal:
    # Whole time points are assigned in chronological order: train < validation < test
    n_tr = min(max(1, round(TRAIN_FRAC * len(times))), len(times) - 2)
    n_va = min(max(1, round(VAL_FRAC * len(times))), len(times) - n_tr - 1)
    part = t_of.map(lambda t: "train" if t < n_tr else "val" if t < n_tr + n_va else "test")
    # A group (subject/replicate set) stays in its earliest partition; later samples are quarantined
    home = part.map(rank).groupby(meta["group"]).transform("min")
    part[part.map(rank) != home] = "quarantined"
else:
    # Static data: group-aware random split (no group is shared between partitions)
    idx, g = np.arange(len(meta)), meta["group"].to_numpy()
    tr, rest = next(GroupShuffleSplit(1, test_size=1 - TRAIN_FRAC, random_state=SEED).split(idx, groups=g))
    test_size = (1 - TRAIN_FRAC - VAL_FRAC) / (1 - TRAIN_FRAC)
    va, te = next(GroupShuffleSplit(1, test_size=test_size, random_state=SEED).split(rest, groups=g[rest]))
    part = pd.Series(index=meta.index, dtype=object)
    part.iloc[tr], part.iloc[rest[va]], part.iloc[rest[te]] = "train", "val", "test"

ids = {k: list(meta.index[part == k]) for k in rank}
if not ids["train"]:
    raise ValueError("Training partition is empty.")


# 3. Preprocessing, selection and network construction: fitted on the given samples only
def fit_prep(sample_ids):
    d = X.loc[sample_ids]
    cols = d.columns[d.notna().any()]
    imp = SimpleImputer(strategy="median").fit(d[cols])
    var = imp.transform(d[cols]).var(axis=0)
    keep = np.argsort(-var)[:TOP_GENES]
    keep = keep[var[keep] > 0]
    scaler = StandardScaler().fit(imp.transform(d[cols])[:, keep])
    return {"cols": cols, "imp": imp, "keep": keep, "scaler": scaler,
            "genes": list(cols[keep]), "fit_ids": list(sample_ids)}


def apply_prep(prep, sample_ids):
    d = X.loc[sample_ids, prep["cols"]]
    return prep["scaler"].transform(prep["imp"].transform(d)[:, prep["keep"]])


def build_network(sample_ids):
    prep = fit_prep(sample_ids)
    G = nx.Graph()
    G.add_nodes_from(prep["genes"])
    if len(sample_ids) >= MIN_SAMPLES and len(prep["genes"]) > 1:
        R = np.corrcoef(apply_prep(prep, sample_ids), rowvar=False)
        a, b = np.triu_indices_from(R, k=1)
        keep = np.abs(R[a, b]) >= MIN_ABS_CORR
        G.add_edges_from((prep["genes"][i], prep["genes"][j], {"correlation": float(R[i, j])})
                         for i, j in zip(a[keep], b[keep]))
    return G, prep


def digest(G):
    edges = sorted((*sorted(e[:2]), round(e[2]["correlation"], 8)) for e in G.edges(data=True))
    return hashlib.sha256(repr((sorted(G.nodes), edges)).encode()).hexdigest()[:16]


# 4. Dynamic updates: one immutable network version per training time point (expanding window)
versions = []
for t in sorted(t_of[ids["train"]].unique()):
    used = [s for s in ids["train"] if t_of[s] <= t]
    G, prep = build_network(used)
    versions.append({"version": len(versions) + 1, "time": times[t], "t": int(t),
                     "samples": used, "graph": G, "prep": prep, "digest": digest(G)})

# 5. Leakage detection
violations = []
h = pd.Series(pd.util.hash_pandas_object(X.round(6), index=False).to_numpy(), index=X.index)
for a, b in [("train", "val"), ("train", "test"), ("val", "test")]:
    A, B = ids[a], ids[b]
    if set(A) & set(B):
        violations.append(f"sample overlap: {a}/{b}")
    if set(meta.loc[A, "group"]) & set(meta.loc[B, "group"]):
        violations.append(f"group/subject overlap: {a}/{b}")
    if set(h[A]) & set(h[B]):
        violations.append(f"identical expression profile in {a} and {b}")
    if temporal and A and B and t_of[A].max() >= t_of[B].min():
        violations.append(f"temporal order violated: {a} is not strictly earlier than {b}")
for v in versions:
    if not set(v["samples"]) <= set(ids["train"]):
        violations.append(f"version {v['version']}: non-training samples used")
    if t_of[v["samples"]].max() > v["t"]:
        violations.append(f"version {v['version']}: future observations used")
    if not set(v["prep"]["fit_ids"]) <= set(v["samples"]):
        violations.append(f"version {v['version']}: preprocessing fitted outside its window")

# 6. Pipeline validation: history is unchanged and rebuilding reproduces every version
immutable = all(digest(v["graph"]) == v["digest"] for v in versions)
reproducible = all(digest(build_network(v["samples"])[0]) == v["digest"] for v in versions)
if not immutable:
    violations.append("historical network version was modified")

# Held-out data are only transformed with training-fitted preprocessing
final = versions[-1]
held_out = {k: apply_prep(final["prep"], ids[k]) for k in ("val", "test") if ids[k]}

# Leakage-prone reference (all samples used for preprocessing and correlations): comparison only
G_leaky, _ = build_network(list(meta.index))
E_train, E_leaky = ({frozenset(e) for e in g.edges} for g in (final["graph"], G_leaky))
jaccard = len(E_train & E_leaky) / len(E_train | E_leaky) if E_train | E_leaky else 1.0

# 7. Audit log
audit = {
    "config": {"train_frac": TRAIN_FRAC, "val_frac": VAL_FRAC, "top_genes": TOP_GENES,
               "min_abs_corr": MIN_ABS_CORR, "min_samples": MIN_SAMPLES, "seed": SEED},
    "partition_mode": "temporal" if temporal else "group-aware random",
    "partitions": ids,
    "quarantined_samples": list(meta.index[part == "quarantined"]),
    "versions": [{"version": v["version"], "time": v["time"], "samples": v["samples"],
                  "n_nodes": v["graph"].number_of_nodes(), "n_edges": v["graph"].number_of_edges(),
                  "n_selected_genes": len(v["prep"]["genes"]), "digest": v["digest"]} for v in versions],
    "violations": violations,
    "history_immutable": immutable,
    "reproducible": reproducible,
    "leaky_reference_edge_jaccard": jaccard,
}
with open(AUDIT_OUTPUT, "w") as f:
    json.dump(audit, f, indent=2, default=str)

print("Partition mode:", audit["partition_mode"])
print({k: len(v) for k, v in ids.items()}, "| quarantined:", len(audit["quarantined_samples"]))
for v in audit["versions"]:
    print(f"Version {v['version']} (time <= {v['time']}): {len(v['samples'])} samples, "
          f"{v['n_nodes']} genes, {v['n_edges']} edges, digest {v['digest']}")
print("Held-out matrices (train-fitted transform):", {k: m.shape for k, m in held_out.items()})
print("Leakage violations:", violations if violations else "none detected")
print("History immutable:", immutable, "| reproducible:", reproducible)
print(f"Edge Jaccard, controlled vs leakage-prone network: {jaccard:.3f}")
print("Audit log written to", AUDIT_OUTPUT)
