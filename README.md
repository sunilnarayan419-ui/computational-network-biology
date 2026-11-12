# Computational Network Biology

### A Collaborative Research and Computational Implementation Initiative

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.12%2B-blue.svg)](https://www.python.org/)
[![Research](https://img.shields.io/badge/Focus-Network%20Biology-purple.svg)]()
[![Open Science](https://img.shields.io/badge/Approach-Open%20Research-green.svg)]()

**Computational Network Biology** is a structured research initiative for understanding biological systems through network science, graph theory, mathematical modelling, and computational analysis.

The project bridges published network biology research and practical computational implementation: it studies the literature, formulates research questions, and translates selected concepts into reproducible computational solutions.

The repository is organized around **100 research-paper workspaces across five interconnected domains**, covering classical network topology, molecular interaction networks, network medicine, higher-order biological networks, and graph-based machine learning.

---

## Table of Contents

1. [At a Glance](#at-a-glance)
2. [Research Objectives](#research-objectives)
3. [Research Domains](#research-domains)
4. [Research Methodology](#research-methodology)
5. [Repository Structure](#repository-structure)
6. [Paper Workspace Standard](#paper-workspace-standard)
7. [Getting Started](#getting-started)
8. [Computational Technologies](#computational-technologies)
9. [Scientific Quality Principles](#scientific-quality-principles)
10. [Collaboration and Contributions](#collaboration-and-contributions)
11. [Current Development Status](#current-development-status)
12. [Research Resources](#research-resources)
13. [Citation and Attribution](#citation-and-attribution)
14. [License](#license)
15. [Community and Acknowledgements](#community-and-acknowledgements)

---

## At a Glance

| | |
|---|---|
| **Focus** | Network biology, graph theory, computational modelling |
| **Scale** | 100 planned paper workspaces in 5 domains |
| **Language** | Python 3.12+ |
| **License** | MIT |
| **Contributions** | Open, branch-based pull-request workflow |
| **Community** | [WhatsApp Research Group](https://chat.whatsapp.com/IbDYrNAFG5FGb6yidDC9yI) |

---

## Research Objectives

- Build a structured understanding of computational network biology through primary scientific literature.
- Translate biological concepts and theoretical frameworks into clearly defined computational problems.
- Explore the mathematical and graph-theoretical foundations of biological network organization.
- Implement selected methods in Python using relevant scientific libraries.
- Provide reproducible implementations with transparent assumptions, limitations, and explanations.
- Encourage collaborative discussion, peer review, and knowledge exchange.
- Identify computational limitations and opportunities for further investigation.

> The emphasis is on scientific understanding and computational reasoning, not on reproducing code from existing publications.

---

## Research Domains

The collection is divided into five domains.

| Domain | Research Area | Paper Range |
|:---:|---|:---:|
| 01 | Classic Foundations & Network Topology | 001–020 |
| 02 | Network Motifs, Modules & Community Detection | 021–040 |
| 03 | Network Medicine & Disease Interactomes | 041–060 |
| 04 | Dynamic, Spatial & Higher-Order Topologies | 061–080 |
| 05 | Graph Representation Learning, Causal Inference & Integration | 081–100 |

<details>
<summary><b>Domain 01: Classic Foundations & Network Topology</b></summary>

- Graph theory and biological network representation
- Degree distributions and network connectivity
- Centrality measures and network hubs
- Scale-free network organization
- Network motifs and structural properties
- Network robustness and perturbation analysis
- Random network models

</details>

<details>
<summary><b>Domain 02: Network Motifs, Modules & Community Detection</b></summary>

- Network motif identification
- Functional modules and molecular complexes
- Community detection algorithms
- Modularity optimization
- Protein interaction networks
- Clustering and network organization
- Functional interpretation of network communities

</details>

<details>
<summary><b>Domain 03: Network Medicine & Disease Interactomes</b></summary>

- Disease-associated molecular networks
- Disease modules and disease-gene prioritization
- Network medicine approaches
- Drug-target interaction networks
- Network propagation algorithms
- Molecular comorbidity
- Network-based therapeutic discovery

</details>

<details>
<summary><b>Domain 04: Dynamic, Spatial & Higher-Order Topologies</b></summary>

- Dynamic biological networks
- Temporal and multilayer networks
- Hypergraphs and higher-order interactions
- Single-cell network analysis
- Spatial biological networks
- Chromatin interaction networks
- Context-dependent network organization

</details>

<details>
<summary><b>Domain 05: Graph Representation Learning, Causal Inference & Integration</b></summary>

- Graph embeddings and representation learning
- Graph neural networks (GNNs)
- Biological network reconstruction
- Causal network inference
- Multi-omics network integration
- Predictive modelling of molecular interactions
- Machine learning for complex biological systems

</details>

---

## Research Methodology

Every paper workspace follows the same seven-stage workflow.

| Stage | Name | Description |
|:---:|---|---|
| 1 | Literature Analysis | Study the publication; identify objectives, theory, methods, and findings. |
| 2 | Conceptual Formulation | Develop questions on the biological interpretation, mathematical reasoning, and concepts. |
| 3 | Computational Problem Formulation | Translate concepts into computational questions linked to the conceptual ones. |
| 4 | Mathematical & Computational Analysis | Choose representations, algorithms, methods, and libraries. |
| 5 | Implementation | Write concise Python implementations. |
| 6 | Scientific Interpretation | Explain approach, assumptions, limitations, and biological relevance. |
| 7 | Collaborative Review | Review work, identify methodological weaknesses, document open questions. |

---

## Repository Structure

```text
computational-network-biology/
├── Domain/
│   ├── README.md
│   ├── 01-classic-foundations-network-topology/
│   │   ├── README.md
│   │   └── paper-001/
│   │       ├── README.md
│   │       ├── CITATION.cff
│   │       ├── conceptual_questions/
│   │       └── computational_questions/
│   │           ├── problem_1/
│   │           │   ├── problem.md
│   │           │   ├── solution.py
│   │           │   └── solution_explaination.md
│   │           └── problem_2/
│   │               ├── problem.md
│   │               ├── solution.py
│   │               └── solution_explaination.md
│   ├── 02-network-motifs-modules-community-detection/
│   ├── 03-network-medicine-disease-interactomes/
│   ├── 04-dynamic-spatial-higher-order-topologies/
│   └── 05-graph-representation-learning-causal-inference-integration/
├── doc/
│   └── assets/
├── scripts/
│   └── validate_structure.py
├── .github/
│   └── workflows/
│       └── ci.yml
├── CONTRIBUTING.md
├── LICENSE
├── README.md
└── stru.py
```

> The structure accommodates the planned collection. Individual paper workspaces are populated progressively.

---

## Paper Workspace Standard

Each research paper has its own directory containing:

| File / Directory | Purpose |
|---|---|
| `README.md` | Research overview, objectives, and collaborative workflow |
| `CITATION.cff` | Citation metadata for the original publication |
| `conceptual_questions/` | Conceptual questions and scientific discussions |
| `computational_questions/` | Computational problems and implementations |

### Computational Problem Structure

Each problem uses a minimal three-file layout:

```text
problem_N/
├── problem.md
├── solution.py
└── solution_explaination.md
```

| File | Contents |
|---|---|
| `problem.md` | Scientific background, computational objective, mathematical formulation, inputs, outputs, constraints |
| `solution.py` | Python implementation using appropriate scientific libraries |
| `solution_explaination.md` | Algorithm, mathematical reasoning, implementation choices, assumptions, limitations |

> The problem statement filename varies in some existing contributions (e.g. `problem_statement.md`). The intended standard is `problem.md`.

---

## Getting Started

### Prerequisites

- Git
- Python 3.12 or a compatible environment
- Basic familiarity with Python and graph theory

### 1. Clone the repository

```bash
git clone https://github.com/sunilnarayan419-ui/computational-network-biology.git
cd computational-network-biology
```

### 2. Explore the domains

```bash
cd Domain
```

Navigate to the domain and paper workspace that matches your interests.

### 3. Run a solution

Each solution is intended to be independently understandable and executable:

```bash
python Domain/01-classic-foundations-network-topology/paper-001/computational_questions/problem_1/solution.py
```

Install any problem-specific dependencies documented in its explanation file before running.

### 4. Validate the repository structure

```bash
python scripts/validate_structure.py
```

A GitHub Actions workflow also runs structural validation and Python syntax compilation.

---

## Computational Technologies

Python is the primary implementation language. Libraries are chosen per problem; use the smallest appropriate set of dependencies.

| Library | Potential Application |
|---|---|
| NumPy | Numerical computing and matrix operations |
| SciPy | Scientific algorithms and statistical analysis |
| Pandas | Biological data manipulation |
| NetworkX | Graph representation and network analysis |
| Matplotlib | Scientific visualization |
| scikit-learn | Statistical learning and machine learning |
| Biopython | Biological sequence analysis |

These are suggestions, not mandatory dependencies.

---

## Scientific Quality Principles

All contributions should follow these principles:

1. **Scientific traceability:** Connect questions and implementations to identifiable publications.
2. **Mathematical clarity:** Define formulations, variables, assumptions, and objectives.
3. **Reproducibility:** Give enough information to understand and reproduce procedures.
4. **Minimal implementation:** Prefer concise, readable code over unnecessary abstraction.
5. **Biological relevance:** Interpret outputs in their biological context.
6. **Critical evaluation:** Identify limitations, uncertainties, and sources of error.
7. **Methodological distinction:** Separate methods from the original publication from independent exercises or extensions.
8. **Responsible research:** Do not treat computational predictions or statistical associations as experimentally established mechanisms.

---

## Collaboration and Contributions

Students, researchers, computational biologists, bioinformaticians, and developers are welcome. You can contribute through:

- Literature analysis and paper selection
- Conceptual question development
- Computational problem formulation
- Mathematical and algorithmic analysis
- Python implementation
- Scientific interpretation and documentation
- Peer review and methodological discussion

### Research Roles

| Role | Responsibilities |
|---|---|
| **Researcher A: Conceptual Research** | Develop 10 conceptual questions; examine theoretical and biological foundations; prepare interpretations and discussion points. |
| **Researcher B: Computational Research** | Develop about 12 computational questions; link them to the conceptual questions; identify mathematical methods and libraries; write concise implementations and explanations. |

Both researchers review the complete work and document disagreements, assumptions, and limitations.

### Contribution Workflow

Please read [CONTRIBUTING.md](CONTRIBUTING.md) first.

1. Create a focused working branch.
2. Work within the relevant paper workspace.
3. Commit with descriptive messages.
4. Push the branch to GitHub.
5. Open a Pull Request for review.
6. Address review comments before merging.

---

## Current Development Status

| Component | Current Scope |
|---|---|
| Research domain framework | 5 domains organized |
| Literature collection | 100 planned paper workspaces |
| Initial paper workspace | Established |
| Computational problems | Initial set of 13 problem directories |
| Conceptual research | Collaborative development |
| Repository validation | GitHub Actions workflow configured |
| Contributions | Open |

> A planned workspace does not imply that its literature review, computational formulation, or implementation is complete.

---

## Research Resources

**Literature sources**

- [Nature](https://www.nature.com/) · [Nature Reviews Genetics](https://www.nature.com/nrg/)
- [PubMed](https://pubmed.ncbi.nlm.nih.gov/) · [Europe PMC](https://europepmc.org/)
- [Google Scholar](https://scholar.google.com/)
- [bioRxiv](https://www.biorxiv.org/) · [arXiv](https://arxiv.org/)

**Computational references**

- [NetworkX](https://networkx.org/documentation/stable/)
- [SciPy](https://docs.scipy.org/doc/scipy/)
- [NumPy](https://numpy.org/doc/)
- [Biopython](https://biopython.org/)

---

## Citation and Attribution

Contributions should cite the original publications they are based on. Where available, each paper workspace includes a `CITATION.cff` file with bibliographic metadata.

Example (initial paper workspace):

> Barabási, A.-L., & Oltvai, Z. N. (2004). Network biology: Understanding the cell's functional organization. *Nature Reviews Genetics*, 5, 101–113. https://doi.org/10.1038/nrg1272

This project does not claim ownership of the original publications. Implementations, explanations, and independent exercises should be clearly distinguished from the original authors' work.

---

## License

Distributed under the [MIT License](LICENSE). It applies to the repository's original code and materials as specified by the license. Original publications and third-party materials remain subject to their own copyright and licensing terms.

---

## Community and Acknowledgements

This initiative is built on collaborative scientific learning, computational experimentation, and open research practices. Contributions, critical feedback, and methodological discussions from the community are welcome.

**Join the research community:** discuss papers, share resources, collaborate on implementations, and review contributions in our WhatsApp group.

[**Join Our WhatsApp Community →**](https://chat.whatsapp.com/IbDYrNAFG5FGb6yidDC9yI)

**Repository:** [Computational Network Biology](https://github.com/sunilnarayan419-ui/computational-network-biology)

*Let's connect biological complexity with computational reasoning, one research problem at a time.*