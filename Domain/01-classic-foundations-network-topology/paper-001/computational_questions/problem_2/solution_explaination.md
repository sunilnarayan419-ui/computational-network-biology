# Solution 02: Data Leakage Prevention in Dynamic Biological Network Pipelines

## 1. Solution Overview

This solution implements a compact, leakage-controlled pipeline for dynamic gene-expression networks. Every step that learns something from data (imputation, gene selection, scaling, correlation, network construction) is executed only on samples that are permitted for that step. Validation and test samples are never used to build or update a network; they are only transformed with preprocessing fitted on training data.

The script follows seven stages:

1. Input loading and validation.
2. Temporal or group-aware partitioning.
3. Preprocessing fitted on permitted samples only.
4. Versioned dynamic network construction.
5. Leakage detection.
6. Pipeline validation (immutability and reproducibility).
7. Audit logging.

## 2. Technology Selection

| Library | Computational responsibility |
|---|---|
| Pandas | Long-to-wide reshaping, sample metadata, partition bookkeeping, row hashing |
| NumPy | Correlation matrix and thresholding |
| scikit-learn | Median imputation, standardization, group-aware splitting |
| NetworkX | Graph representation of each network version |

`hashlib`, `json` and `re` from the Python standard library are used for network fingerprints, the audit file and natural sorting of time labels. No MLOps framework is used.

## 3. Computational Methodology

### Stage 1: Input Loading and Validation

**Library:** Pandas

The input is a long-format table with columns `gene`, `sample`, `time`, `expression` and an optional `group` column (subject, donor or replicate set). If `group` is absent, each sample is treated as its own group.

The script stops with an error if required columns are missing, if `gene`, `sample` or `time` contain missing values, if a gene-sample pair is duplicated, or if one sample carries more than one time point or group. The table is pivoted to a samples × genes matrix `X`.

### Stage 2: Dataset Partitioning

**Library:** Pandas, scikit-learn

Two modes are used depending on the data:

- **Temporal mode (three or more distinct time points).** Time points are ordered naturally (`T2` before `T10`) and whole time points are assigned in chronological order to training, validation and test. The number of time points per partition follows `TRAIN_FRAC` and `VAL_FRAC`, with at least one time point each. Because a group can span several time points, each group is kept in its earliest partition. Its later samples are labelled `quarantined` and excluded from all partitions.
- **Static mode (fewer than three time points).** `GroupShuffleSplit` creates a random split in which no group appears in more than one partition. The split is reproducible through `SEED`.

### Stage 3: Preprocessing Isolation

**Library:** scikit-learn, NumPy

The function `fit_prep(sample_ids)` fits all learned transformations using only the supplied samples:

1. Genes with no measured value are dropped.
2. A median imputer is fitted.
3. The `TOP_GENES` genes with the largest variance (and variance > 0) are selected.
4. A standard scaler is fitted on the selected genes.

`apply_prep()` applies the fitted objects to any sample set without re-fitting. Validation and test matrices are produced only by `apply_prep()` with the final training-fitted preprocessing.

### Stage 4: Versioned Dynamic Network Construction

**Library:** NumPy, NetworkX

For each training time point `t`, a network version is built from all training samples with time rank ≤ `t` (an expanding window). Each version has its own preprocessing fit, so a historical state never contains information from later time points.

For the selected genes, the Pearson correlation matrix is computed on the standardized training samples:

r(i,j) = Cov(Zi, Zj) / (σi × σj)

An undirected edge is created when |r(i,j)| ≥ `MIN_ABS_CORR`. Networks are built only when at least `MIN_SAMPLES` samples are available; otherwise a version contains nodes and no edges. The threshold is a fixed configuration value and is never selected using validation or test data.

Each version stores its sample list, fitted preprocessing, graph and a SHA-256 fingerprint of its nodes and rounded edge weights. New versions are added; existing versions are not modified.

### Stage 5: Leakage Detection

**Library:** Pandas

For every pair of partitions the script checks:

| Check | Leakage category |
|---|---|
| Shared sample IDs | Sample leakage |
| Shared groups | Sample / subject leakage |
| Identical rounded expression profiles in different partitions | Duplicated records |
| Latest time in the earlier partition not strictly before the earliest time in the later one (temporal mode only) | Temporal leakage |

For every network version it checks that all samples come from the training partition, that no sample is later than the version's time point, and that the preprocessing fit used only the version's own samples (preprocessing, association and network leakage).

Evaluation leakage is prevented by design: no step of the script reads validation or test values to choose thresholds, gene counts or parameters. This property is not detected at run time.

### Stage 6: Pipeline Validation

**Library:** NetworkX, standard library

- **Immutability:** the fingerprint of each stored graph is recomputed and compared with the fingerprint recorded at creation.
- **Reproducibility:** each version is rebuilt from its recorded sample list and its fingerprint is compared with the stored one.
- **Leakage-prone reference:** a network is also built using all samples (including validation, test and quarantined samples). The Jaccard index between its edge set and the final controlled network quantifies how much the leaky pipeline changes the network. This reference is for comparison only.

### Stage 7: Audit Log

**Library:** standard library

`leakage_audit.json` records configuration, partition mode, sample lists per partition, quarantined samples, every network version (time point, samples, node and edge counts, number of selected genes, fingerprint), detected violations, the immutability and reproducibility results, and the Jaccard comparison.

## 4. Important Mathematical Operations

**Median imputation:** a missing value of gene g is replaced by the median of g in the fitting samples.

**Standardization:** Z = (X − μ_train) / σ_train, where μ and σ come from the fitting samples only.

**Variance selection:** genes are ranked by variance of the imputed fitting data.

**Pearson correlation:** as given in Stage 4. Because Pearson correlation is invariant to centering and scaling, standardization does not change the correlations. It is retained because the same fitted scaler is required for held-out data transformation.

**Edge Jaccard index:** J = |E_controlled ∩ E_leaky| / |E_controlled ∪ E_leaky|.

## 5. Input Data Requirements

| gene | sample | time | expression | group (optional) |
|---|---|---|---|---|
| Gene A | S01 | T1 | 2.45 | Donor1 |
| Gene B | S01 | T1 | 3.12 | Donor1 |
| Gene A | S02 | T2 | 2.87 | Donor1 |

- One row per gene-sample measurement; no duplicate gene-sample pairs.
- One time point and one group per sample.
- Numeric expression values; non-numeric entries are treated as missing.
- Time labels sort naturally (`T1`, `T2`, `T10`, or `1`, `2`, `10`).

No data file is supplied with the repository; the user provides `expression_long.csv`.

## 6. Output Interpretation

| Output | Meaning |
|---|---|
| Partition sizes | Samples in training, validation and test; quarantined count |
| Network versions | Nodes, edges and fingerprint for each time window |
| Held-out matrix shapes | Validation and test data transformed with training-fitted preprocessing |
| Leakage violations | List of detected problems; "none detected" means none of the implemented checks fired |
| History immutable / reproducible | Results of the Stage 6 checks |
| Edge Jaccard | Similarity between the controlled and the leakage-prone network; a low value shows that leakage would change the network |
| `leakage_audit.json` | Complete record of partitions and network history |

## 7. Assumptions and Limitations

1. The checks cover the leakage routes listed in Stage 5. Leakage that comes from outside the script (for example, genes pre-filtered on the whole dataset before the file was created) cannot be detected.
2. Duplicate detection compares expression profiles rounded to six decimals. Near-duplicates and technical replicates without a shared `group` value are not detected.
3. Temporal mode assigns whole time points to partitions. With few time points, partitions can be very small, and early versions may contain no edges.
4. In static mode, the time-point windows used for versions are not chronological evidence, because the split is random.
5. Quarantining the later samples of a group reduces the data available for validation and test.
6. Correlation edges are statistical associations, not confirmed regulatory or physical interactions. The script applies a correlation threshold only, without significance testing or multiple-testing correction.
7. The median imputer and variance selection are simple choices; other preprocessing can be added only if it is fitted inside `fit_prep()`.
8. The script does not train a predictive model, so "performance inflation" is shown only through the network-difference (Jaccard) comparison.
