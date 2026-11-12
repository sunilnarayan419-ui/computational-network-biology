# Solution 09: Identification of Highly Connected Nodes and Preferential Attachment

## 1. Solution Overview

This solution does two separate jobs. The first identifies and ranks the most highly connected nodes in a network and highlights candidate hubs. The second, which runs only when time-stamped edges are supplied, tests whether new connections are directed preferentially toward nodes that already have many connections.

The two parts are deliberately kept apart. Finding a hub in a final network does not show that preferential attachment generated it. Without time information the script stops after hub identification and topology description.

The implementation follows five computational stages:

1. Edge-list loading and cleaning.
2. Degree calculation, maximum-degree set, ranking, and hubs.
3. Other centrality measures, reported separately.
4. Degree-distribution and hub visualization.
5. Maximum-likelihood estimation of the attachment exponent from time-stamped edges.

## 2. Technology Selection

| Library | Computational responsibility |
|---|---|
| Pandas | Edge list, ranked tables, annotations |
| NumPy | Degree histograms and likelihood calculations |
| NetworkX | Graph construction, degree, betweenness, eigenvector centrality |
| SciPy | Spearman correlation, likelihood optimisation, root finding, chi-square test |
| Matplotlib | Degree distribution and hub-highlighted network |

Statsmodels is not used. The attachment exponent is estimated by direct maximum likelihood with SciPy, which handles the changing pool of eligible nodes correctly, so a regression library is not required.

## 3. Computational Methodology

### Stage 1: Loading and Cleaning

**Library:** Pandas

The input is an edge list with `source` and `target` columns and an optional numeric `time` column. Self-loops and missing records are removed. For undirected networks, A–B and B–A are the same edge. Duplicate edges are merged and the earliest occurrence is kept. Nodes listed in `EXCLUDE_NODES` (for example currency metabolites) are removed, which allows a filtering sensitivity check.

For metabolic data, the node type must be fixed by the input. A reaction–metabolite bipartite edge list gives each metabolite a degree equal to the number of reactions it takes part in. Projections can inflate the degree of currency metabolites.

### Stage 2: Degree, Maximum-Degree Set, and Hubs

**Library:** NetworkX, Pandas

For undirected networks the degree $k_i$ is the number of incident edges. For directed networks, in-degree $k_{\text{in}}(i)$ and out-degree $k_{\text{out}}(i)$ are calculated and handled as two separate measures. They are never added together.

The maximum degree and the highest-connectivity set are:

$$k_{\max} = \max_i k_i, \qquad H = \{\, i \in V : k_i = k_{\max} \,\}$$

All nodes tied at $k_{\max}$ are listed. Nodes are ranked in descending degree with tied nodes sharing a rank. Candidate hubs are all nodes whose degree is at least the `HUB_QUANTILE` quantile (default 0.95), ties included. This threshold is a screening choice and not a biological definition of a hub.

### Stage 3: Other Centralities

**Library:** NetworkX, SciPy

Betweenness centrality (sampled when the network exceeds `BETWEENNESS_SAMPLES` nodes) and eigenvector centrality are added as separate columns. The Spearman correlation between degree and betweenness is printed to show that degree and betweenness are related but not interchangeable. Eigenvector centrality may be undefined for some directed or disconnected networks, and the script reports when it could not be computed. An optional annotation file (`node` column plus any annotation columns) is merged into the hub table.

### Stage 4: Visualization

**Library:** Matplotlib

The figure `hubs.png` shows each degree distribution on log–log axes with the hub threshold marked. Networks with at most `MAX_DRAW_NODES` nodes are also drawn with candidate hubs in red. The degree distribution is saved to `degree_distribution_<measure>.csv` and the full ranked table to `hub_table.csv`.

### Stage 5: Preferential-Attachment Test

**Library:** NumPy, SciPy

This stage requires a `time` column. Edges are processed in time order, and edges sharing a time stamp are treated as one snapshot.

An attachment event is an edge that joins a node new to the network to a node already present. For directed networks the new node must be the source, and the chosen node is the target, so the measure is the in-degree of the target. For each event the degree of the chosen node is recorded before that time step is applied. Edges between two existing nodes or two new nodes are not counted as events.

The model for the probability that a given event selects existing node $i$ is:

$$\Pi(i) = \frac{(k_i + a)^{\alpha}}{\sum_j (k_j + a)^{\alpha}}$$

where the sum covers all nodes existing before the event, $a \ge 0$ is the initial-attractiveness offset (`INITIAL_ATTRACTIVENESS`), and $\alpha$ is the attachment exponent. With $a = 0$, nodes of degree 0 cannot be selected and events that select such nodes are not counted.

The log-likelihood over all events is:

$$\ell(\alpha) = \alpha \sum_{\text{events}} \ln(k_{\text{chosen}} + a) - \sum_{\text{events}} \ln \sum_j (k_j + a)^{\alpha}$$

The denominator uses the degree histogram of existing nodes at each time step, so the changing number and degrees of eligible nodes are accounted for. $\alpha$ is found by bounded numerical maximisation of $\ell$ between $-1$ and $4$. A 95% interval is taken from the profile likelihood, using the points where $\ell$ falls by 1.92 below its maximum.

Interpretation of $\alpha$:

- $\alpha \approx 1$: linear preferential attachment.
- $\alpha > 1$: superlinear attachment.
- $0 < \alpha < 1$: sublinear attachment.
- $\alpha \approx 0$: attachment independent of degree.

Likelihood-ratio tests compare the fitted model with uniform attachment ($\alpha = 0$) and linear preferential attachment ($\alpha = 1$). The statistic $2(\ell_{\max} - \ell(\alpha_0))$ is compared with a chi-square distribution with one degree of freedom. Observed attachment counts per degree class are compared with counts expected under the uniform, linear, and fitted models and saved to `attachment_pattern.csv`.

If fewer than 30 events are available, no exponent is estimated.

## 4. Input Data Requirements

| Item | Requirement |
|---|---|
| `INPUT` | CSV with `source`, `target`, and optionally a numeric `time` column |
| `DIRECTED` | Edge direction convention (new node as source for attachment testing) |
| `ANNOTATIONS` | Optional CSV with a `node` column |
| `EXCLUDE_NODES` | Optional nodes to remove for filtering sensitivity |

Time stamps may be real dates converted to numbers, sequence indices, or snapshot numbers. Edge order must reflect the order in which connections were acquired, which is rarely certain for biological data.

## 5. Expected Computational Output

| Output | Description |
|---|---|
| Node and edge counts | Size of the cleaned network |
| $k_{\max}$ and tied nodes | Highest-connectivity set for each degree measure |
| Ranked table and candidate hubs | `hub_table.csv` with centralities and annotations |
| Degree distribution | `degree_distribution_<measure>.csv` and `hubs.png` |
| $\alpha$, confidence interval, likelihood-ratio tests | Attachment-exponent results when time data are present |
| Observed vs expected attachment | `attachment_pattern.csv` |

Without a time column, the script prints that preferential attachment cannot be tested and limits its conclusions to hub identification.

## 6. Scientific Assumptions and Limitations

1. A high-degree node does not show that preferential attachment generated the network.
2. Estimates assume that file order or time stamps reflect real growth order, and that nodes enter the data when they first appear in an edge.
3. Missing, delayed, or database-biased observations bias the estimate. Changes in database coverage over time are not corrected.
4. Only a power-law kernel $(k + a)^\alpha$ is fitted. Other kernels and alternative growth models such as fitness or duplication models are not compared beyond the uniform and linear cases.
5. Events within one time stamp are treated as sharing the same pre-snapshot degrees, which is an approximation for large snapshots.
6. Weights and node types are not used, and the alternative growth mechanisms listed in the original problem are not simulated.
7. A hub is not the same as a node that drives network growth.
8. Snapshot histograms are stored for every time step with events, so memory use grows with the number of time steps and the maximum degree.

## 7. Conclusion

The script identifies and ranks highly connected nodes, including ties and directed degrees, and, when time-stamped data are available, estimates an attachment exponent with uncertainty and compares it with uniform and linear models. Hub identification and growth-mechanism testing are reported as separate results.
