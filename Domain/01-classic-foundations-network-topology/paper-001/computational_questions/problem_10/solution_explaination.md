# Solution 10: Conditional-Probability-Guided Search for Larger Network Motifs

## 1. Solution Overview

This solution estimates how often a chosen larger motif (3 to 6 nodes) occurs in a biological network without enumerating every node subset. Node sets are grown one node at a time, and conditional probabilities learned from a pilot sample steer the growth toward configurations that complete the target motif.

Because guided sampling selects some configurations more often than others, each sampled configuration is weighted by the inverse of its exact sampling probability. The result is an unbiased frequency estimate with bootstrap uncertainty, not a naive sample proportion.

The script keeps three tasks separate:

1. **Candidate discovery:** distinct node sets matching the motif that the guided search finds.
2. **Frequency estimation:** an importance-sampling estimate of the number and fraction of matching node sets.
3. **Significance testing:** comparison with degree-preserving randomized networks processed by the same pipeline.

The implementation follows six computational stages:

1. Network loading and motif definition.
2. Uniform pilot sampling.
3. Conditional-probability estimation.
4. Probability-guided sequential sampling.
5. Importance-weighted frequency estimate with uncertainty.
6. Null-model comparison with multiple-testing correction.

## 2. Technology Selection

| Library | Computational responsibility |
|---|---|
| Pandas | Edge list, probability tables, output files |
| NumPy | Random sampling, weights, bootstrap |
| NetworkX | Graph construction, induced subgraphs, isomorphism, degree-preserving randomization |
| SciPy | Beta-distribution intervals for conditional probabilities |
| statsmodels | Benjamini–Hochberg correction across motifs |
| Matplotlib | Observed estimate vs null distribution |

## 3. Computational Methodology

### Stage 1: Network and Motif Definition

**Library:** Pandas, NetworkX

The input is an edge list. Self-loops are removed and duplicates are merged. Each target motif is defined in `MOTIFS` as a list of edges on node names of the user's choice. The motif must be weakly connected and have 3 to 6 nodes, and all motifs in one run must have the same size $m$.

A motif $M = (V_M, E_M)$ is matched as an **induced** pattern. A sampled node set $X$ is a match when the subgraph induced on $X$ is isomorphic to $M$, which means exactly the motif's edges are present and no others. Edge direction is preserved for directed networks. Undirected networks are internally represented with an arc in each direction.

Node types, edge weights, and interaction labels are not used. For metabolic data, the user decides whether nodes are metabolites or reactions through the edge list supplied.

### Stage 2: Sequential Construction of Node Sets

**Library:** NetworkX, NumPy

A candidate configuration is a connected set of $m$ nodes built in $m$ steps. The first node is chosen uniformly among nodes with at least one edge. Each later step adds one node from the **frontier** of the current set $S$, which is the set of nodes outside $S$ with at least one edge to or from $S$.

Each frontier node $u$ is described by its extension context:

$$B = (|S|,\; \text{arcs inside } S,\; \text{arcs } u \to S,\; \text{arcs } S \to u)$$

The context summarises the partial pattern and how the candidate attaches to it. It is a coarse description and does not distinguish all non-isomorphic partial patterns.

### Stage 3: Conditional-Probability Estimation

**Library:** NumPy, SciPy

In the pilot, `PILOT_SAMPLES` node sets are built with uniform choices at every step. For each context $B$ visited, the script counts the number of times it occurred, $n_B$, and the number of times the final set was the target motif, $n_{AB}$. The conditional probability is:

$$P(A \mid B) = \frac{\text{count}(A \text{ and } B)}{\text{count}(B)}$$

where $A$ means the finished set matches the motif. The estimate used for guidance is smoothed with a documented Beta prior centred on the overall pilot hit rate $p_0$:

$$\hat P(A \mid B) = \frac{n_{AB} + \kappa\, p_0}{n_B + \kappa}$$

where $\kappa$ is `SMOOTHING`. Contexts never seen in the pilot take the value $p_0$ and are not treated as impossible. For each observed context the script also reports $n_B$, $n_{AB}$, the raw proportion, the smoothed estimate, and a 95% Jeffreys (Beta) interval.

The probability of a whole extension sequence follows the chain rule without assuming independence, because each step is conditioned on the full context of the partial set built so far:

$$P(E_1, \dots, E_q) = P(E_1)\, P(E_2 \mid E_1) \cdots P(E_q \mid E_1, \dots, E_{q-1})$$

If the pilot finds no motif at all, the guided stage falls back to uniform sampling and says so.

### Stage 4: Guided Sampling

**Library:** NumPy

In the guided stage, each step selects frontier node $u$ with probability proportional to:

$$w(u \mid S) = \hat P(A \mid B_u) + \epsilon$$

where $\epsilon$ is `WEIGHT_FLOOR`. The floor keeps every connected configuration reachable, which is needed for an unbiased estimate. The pilot and guided samples are independent draws, and the guidance weights are fixed before guided sampling begins. This is a sequential importance-sampling design. No beam search is used.

### Stage 5: Importance-Weighted Estimation

**Library:** NumPy

A node set $X$ can be built in several orders. Its total sampling probability $p(X)$ is calculated exactly by summing over all construction orders using dynamic programming over subsets of $X$:

$$p(X) = \sum_{\text{orders } \pi} \frac{1}{N_{\text{seed}}} \prod_{j=1}^{m-1} \frac{w(\pi_j \mid S_{j-1})}{\sum_{u \in \text{frontier}(S_{j-1})} w(u \mid S_{j-1})}$$

where $N_{\text{seed}}$ is the number of nodes that can start a set. With $n$ guided samples, the unbiased estimates are:

$$\hat N_{\text{motif}} = \frac{1}{n}\sum_{s=1}^{n} \frac{\mathbf{1}[X_s \text{ matches motif}]}{p(X_s)}, \qquad \hat N_{\text{total}} = \frac{1}{n}\sum_{s=1}^{n} \frac{1}{p(X_s)}$$

Failed draws, which occur when the seed lies in a component smaller than $m$, count as zero in both sums. $\hat N_{\text{motif}}$ estimates the number of matching node sets. $\hat N_{\text{total}}$ estimates the number of connected $m$-node sets. Their ratio $\hat N_{\text{motif}} / \hat N_{\text{total}}$ is the estimated motif frequency among connected $m$-node sets.

Uncertainty comes from a bootstrap that resamples the guided samples. The script reports 95% percentile intervals for the count and the frequency, and an effective sample size $(\sum w)^2 / \sum w^2$ that shows whether a few very large weights dominate. Runtime of the observed-network search is printed to show the cost of the method.

### Stage 6: Null-Model Significance

**Library:** NetworkX, SciPy, statsmodels

Degree-preserving edge swaps create `N_RANDOM` randomized networks (in-degree and out-degree are preserved for directed networks). Pilot, guided sampling, and estimation are repeated on each randomized network with the same settings. For each motif the empirical p-value for enrichment is:

$$p = \frac{1 + \#\{\hat N_{\text{null}} \ge \hat N_{\text{obs}}\}}{1 + N_{\text{random}}}$$

When several motifs are tested, p-values are adjusted by the Benjamini–Hochberg procedure. The null mean and standard deviation are printed. The figure `motif_null_comparison.png` shows the null distribution with the observed estimate and its interval.

## 4. Input Data Requirements

| Item | Requirement |
|---|---|
| `INPUT` | CSV edge list with `source` and `target` columns |
| `DIRECTED` | Directed or undirected interpretation |
| `MOTIFS` | Dictionary of motif name to edge list, same size $m$ for all motifs |
| `PILOT_SAMPLES`, `GUIDED_SAMPLES`, `N_RANDOM`, `SEED` | Sampling budgets and random seed |

## 5. Expected Computational Output

| Output | Description |
|---|---|
| Pilot summary and `conditional_probabilities_<motif>.csv` | Contexts with $n_B$, $n_{AB}$, probabilities, and intervals |
| `candidate_motifs_<motif>.csv` | Distinct matching node sets found by the guided search |
| Estimated motif count and frequency | With bootstrap intervals and effective sample size |
| Null mean, standard deviation, empirical p, BH-adjusted p | Statistical comparison with randomized networks |
| Runtime | Time for the observed search |
| `motif_null_comparison.png` | Observed vs null estimate |

A motif that is sampled often is not significant unless its estimated count exceeds what the randomized networks give. A small empirical p-value is limited by the number of null networks: with 20 networks the smallest possible value is about 0.048.

## 6. Scientific Assumptions and Limitations

1. Conditional probabilities reduce sampling effort but do not remove combinatorial complexity. The pilot, exact probability calculation, and null networks all have costs, and the probability calculation grows with the number of construction orders ($2^m$ subsets).
2. The conditional-probability intervals treat overlapping pilot samples as independent, so they are somewhat optimistic.
3. The context used for conditioning is coarse, and a better context would improve guidance but not change unbiasedness.
4. The estimator is unbiased given its sampling design, but variance can be large when a few configurations receive very small sampling probabilities. Check the effective sample size.
5. The bootstrap treats guided samples as independent draws. The null estimates themselves contain sampling noise.
6. Degree-preserving randomization is one possible null model. It does not preserve other properties such as clustering or node types.
7. Edge weights, labels, and node types are ignored, and motifs larger than six nodes are not supported.
8. The original validation strategy (comparison with exact enumeration on small networks, coverage checks, uniform vs guided sampling comparisons) is not part of the script and would have to be run separately.

## 7. Conclusion

The script uses a pilot sample to estimate conditional probabilities of motif extensions, guides sequential sampling with them, and corrects for the resulting unequal sampling probabilities with exact inverse-probability weights. Estimated motif frequencies are accompanied by bootstrap uncertainty and compared with degree-preserving randomized networks. Discovery, frequency estimation, and significance are reported as three separate results.
