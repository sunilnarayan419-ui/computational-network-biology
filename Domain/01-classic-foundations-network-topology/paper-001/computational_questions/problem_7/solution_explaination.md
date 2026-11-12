# Solution 07: Statistical Identification of Scale-Free and Non-Scale-Free Networks

## 1. Solution Overview

This solution tests whether the degree distribution of a biological network is statistically consistent with a power law, and compares the power law with competing distributions before reaching a conclusion.

The script never classifies a network from the appearance of a log–log plot. The power-law exponent and lower cutoff are estimated by maximum likelihood, goodness of fit is assessed with a bootstrap Kolmogorov–Smirnov test, and the power law is compared with exponential, lognormal, and truncated power-law models using likelihood ratios.

The implementation follows eight computational stages:

1. Edge-list loading and graph construction.
2. Degree calculation.
3. Maximum-likelihood power-law fitting and $k_{\min}$ estimation.
4. Bootstrap goodness-of-fit testing.
5. Alternative-distribution comparison.
6. Evidence-based classification.
7. Statistical report.
8. Degree-distribution and CCDF visualization.

## 2. Technology Selection

| Library | Computational responsibility |
|---|---|
| Pandas | Loading the edge list and building the degree-frequency table |
| NumPy | Degree arrays, random numbers, and synthetic bootstrap data |
| NetworkX | Graph construction and degree calculation |
| powerlaw | Maximum-likelihood fitting, $k_{\min}$ selection, and likelihood-ratio comparisons |
| Matplotlib | Degree-distribution and CCDF figures |

SciPy is not imported directly because the `powerlaw` package already uses it internally for the numerical work. Weighted-network and optimisation libraries are excluded.

## 3. Computational Methodology

### Stage 1: Graph Construction

**Library:** Pandas, NetworkX

The input is an edge list with `source` and `target` columns. Missing values and self-loops are removed, and duplicate edges collapse into one edge. Undirected networks produce one degree per node. Directed networks (`DIRECTED = True`) are analysed twice, once for in-degree and once for out-degree, and are never merged.

### Stage 2: Degree Sequence

**Library:** NetworkX, NumPy

The degree of node $i$ is:

$$k_i = \sum_j A_{ij}$$

where $A$ is the adjacency matrix. Nodes with $k = 0$ are discarded from the fitted sequence because a discrete power law that starts at a positive degree cannot include them. The degree-frequency table ($N_k$ and $P(k) = N_k / N$) is saved as `degree_frequency_<measure>.csv`.

### Stage 3: Maximum-Likelihood Power-Law Fit

**Library:** powerlaw

The candidate model for degrees $k \ge k_{\min}$ is:

$$P(k) = C k^{-\gamma}$$

For discrete data the exponent $\gamma$ is estimated by maximising the likelihood of the observed tail, not by regression on a log-transformed histogram. The lower cutoff $k_{\min}$ is chosen as the value that minimises the Kolmogorov–Smirnov distance between the empirical and fitted tail distributions:

$$D = \max_{k \ge k_{\min}} \left| F_{\text{empirical}}(k) - F_{\text{fitted}}(k) \right|$$

The script reports $\gamma$, its standard error, $k_{\min}$, the number of nodes in the tail, and $D$.

### Stage 4: Bootstrap Goodness of Fit

**Library:** powerlaw, NumPy

A semi-parametric bootstrap is used. Each synthetic dataset has the same size as the real one. For every synthetic node, with probability $n_{\text{tail}}/N$ the degree is drawn from the fitted power law, and otherwise it is resampled from the observed degrees below $k_{\min}$. Each synthetic dataset is refitted, including a new $k_{\min}$ search, and its KS statistic is recorded.

The bootstrap p-value is the fraction of synthetic datasets whose KS statistic is at least as large as the observed $D$. A small p-value means the power law is rejected. A large p-value means only that the power law has not been rejected.

### Stage 5: Alternative Distributions

**Library:** powerlaw

The power law is compared with these alternatives on the same tail data:

$$\text{Exponential: } P(k) \propto e^{-\lambda k}$$

$$\text{Lognormal: } P(k) \propto \frac{1}{k\sigma\sqrt{2\pi}} \exp\left[-\frac{(\ln k - \mu)^2}{2\sigma^2}\right]$$

$$\text{Truncated power law: } P(k) \propto k^{-\gamma} e^{-\lambda k}$$

For each alternative the script computes the log-likelihood ratio $R$ and its p-value. Positive $R$ favours the power law and negative $R$ favours the alternative. The truncated power law contains the pure power law as a special case, so it is compared with a nested test. The other alternatives use the normalized (Vuong) ratio test.

### Stage 6: Classification

The classification rules are applied in this order:

| Outcome | Condition |
|---|---|
| Inconclusive | Fewer than `MIN_NODES` nodes, or fewer than `MIN_TAIL` observations in the fitted tail, or the bootstrap failed |
| Non-scale-free evidence | Bootstrap p-value below `SIGNIFICANCE`, or an alternative is significantly preferred ($R < 0$ and $p <$ `LR_ALPHA`) |
| Inconclusive | The power law is not rejected, but it is not significantly better than the exponential model either, so a light-tailed explanation cannot be excluded |
| Scale-free supported | The power law is not rejected, the tail and network are large enough, no alternative is significantly preferred, and the exponential model is significantly worse |

The thresholds are configurable screening criteria, not universal constants. Indistinguishable lognormal or truncated power-law fits are reported in the printed comparisons and should be read as a caveat on a "supported" result.

### Stage 7: Visualization

**Library:** Matplotlib

The left panel shows $P(k)$ on log–log axes. The right panel shows the complementary cumulative distribution function (CCDF) with the fitted power law and alternatives over the fitted tail. These plots are descriptive only and are never used to make the classification. The figure is saved as `degree_distribution.png`.

## 4. Input Data Requirements

| Item | Requirement |
|---|---|
| `INPUT` | CSV edge list with columns `source` and `target` |
| `DIRECTED` | `True` for directed networks |
| `BOOTSTRAPS`, `SEED` | Bootstrap replicate count and random seed for reproducibility |
| `MIN_NODES`, `MIN_TAIL` | Size limits below which no classification is made |

An adjacency matrix must be converted to an edge list before use. Edge weights are ignored, so the analysis uses ordinary degree and never node strength. No data are generated by the script.

## 5. Expected Computational Output

| Output | Description |
|---|---|
| Node and edge counts | Size of the cleaned network |
| `degree_frequency_<measure>.csv` | Degree-frequency table |
| $\gamma$, $k_{\min}$, tail size | Fitted power-law parameters |
| $D$ and bootstrap p-value | Goodness-of-fit evidence |
| $R$, p-value, and parameters per alternative | Likelihood-ratio comparisons |
| Classification | Scale-free supported, non-scale-free evidence, or inconclusive |
| `degree_distribution.png` | Log–log degree distribution and CCDF |

Software versions and the random seed are printed so the run can be reproduced. A bootstrap p-value below 0.10 is read as rejection of the power law, following the criterion commonly used with this test.

## 6. Scientific Assumptions and Limitations

1. Passing the goodness-of-fit test does not prove scale-free topology. It only means the power law has not been rejected.
2. The semi-parametric bootstrap assumes that degrees are independent draws, but degrees in a network are not strictly independent.
3. The bootstrap is stochastic, so p-values change slightly with the number of replicates and the seed.
4. Small or heavily sampled networks give little statistical power, which is why they are classified as inconclusive.
5. Database and detection bias can change the apparent degree distribution.
6. Metabolite, reaction, protein, and gene networks differ in construction and should be analysed separately.
7. Isolated nodes are not present in an edge list and degree-zero nodes are excluded from the fit.
8. The script does not generate synthetic validation networks, perform robustness resampling, or fit weighted-network models. Those checks from the original validation strategy would have to be run separately.

## 7. Conclusion

The script combines maximum-likelihood power-law estimation, bootstrap goodness-of-fit testing, and likelihood-ratio comparison with alternatives to classify a network as scale-free supported, non-scale-free evidence, or inconclusive. It is meant to support the scale-free assumptions of earlier problems with statistical evidence, not to replace biological interpretation.
