
"""Create the Computational Network Biology repository structure.

Run from the root of your cloned repository:
    python create_network_biology_structure.py

Safe to rerun: existing files are preserved unless --overwrite is supplied.
"""
from pathlib import Path
import argparse

DOMAINS = [
    (1, "01-classic-foundations-network-topology", "Classic Foundations & Network Topology", 1, 20,
     "Graph theory, centrality, topology, degree distributions, network properties, and structural robustness."),
    (2, "02-network-motifs-modules-community-detection", "Network Motifs, Modules & Community Detection", 21, 40,
     "Motif discovery, clustering, modularity, community detection, and protein-complex identification."),
    (3, "03-network-medicine-disease-interactomes", "Network Medicine & Disease Interactomes", 41, 60,
     "Network propagation, disease modules, disease-gene prioritisation, drug-target networks, and network epidemiology."),
    (4, "04-dynamic-spatial-higher-order-topologies", "Dynamic, Spatial & Higher-Order Topologies", 61, 80,
     "Multilayer networks, hypergraphs, higher-order interactions, single-cell, chromatin, and spatial networks."),
    (5, "05-graph-representation-learning-causal-inference-integration", "Graph Representation Learning, Causal Inference & Integration", 81, 100,
     "Graph embeddings, graph neural networks (GNNs), causal networks, network reconstruction, and multi-omics integration."),
]

ROOT_README = """# Computational Network Biology

A collaborative, paper-driven learning and implementation project by two contributors.

## Scope
We will cover **100 research papers across 5 computational domains**. For every paper:
- Both contributors read the paper and understand its computational context.
- Researcher A prepares **10 conceptual questions** and evidence-based answers.
- Researcher B prepares **12 computational/coding questions** linked to those conceptual questions.
- Both contributors discuss the questions, challenge assumptions, and review answers.
- Relevant methods are implemented, tested, and documented where feasible.

The objective is to understand, critically evaluate, reproduce, and where appropriate extend computational methods used in network biology.

## Domain plan

| Domain | Papers | Count | Computational focus |
|---|---:|---:|---|
| 1. Classic Foundations & Network Topology | 1–20 | 20 | Graph theory, centrality, topology, network properties |
| 2. Network Motifs, Modules & Community Detection | 21–40 | 20 | Motif discovery, clustering, modularity, community detection |
| 3. Network Medicine & Disease Interactomes | 41–60 | 20 | Network propagation, disease modules, drug-target networks |
| 4. Dynamic, Spatial & Higher-Order Topologies | 61–80 | 20 | Multilayer networks, hypergraphs, single-cell and spatial networks |
| 5. Graph Representation Learning, Causal Inference & Integration | 81–100 | 20 | Graph embeddings, GNNs, causal networks, multi-omics integration |
| **Total** | **1–100** | **100** | **5 domains** |

## Weekly workflow
1. Select the next paper from the literature catalogue.
2. Record bibliographic details and DOI/URL.
3. Prepare 10 conceptual questions.
4. Prepare 12 computational/coding questions that build on the conceptual questions.
5. Discuss answers, assumptions, limitations, and alternative approaches.
6. Implement relevant methods and add tests or reproducible examples.
7. Record conclusions and unresolved questions.

## Collaboration
Use a separate branch for each paper or focused contribution. Open a Pull Request, request review from the other contributor, resolve comments, and merge only after repository checks pass. Avoid committing directly to `main`.

See [CONTRIBUTING.md](CONTRIBUTING.md) for the collaboration workflow.
"""

CONTRIBUTING = """# Contributing

## Branching and Pull Requests
1. Keep `main` stable; do not work directly on it.
2. Create a focused branch, e.g. `paper-001/conceptual-questions` or `paper-001/implementation`.
3. Commit with clear, specific messages and push the branch.
4. Open a Pull Request and request review from the other contributor.
5. Resolve comments and ensure CI passes before merging.

## Research responsibilities
- Researcher A drafts 10 conceptual questions.
- Researcher B drafts 12 computational/coding questions explicitly linked to the conceptual questions.
- Both contributors review the complete set and discuss answers.
- Record disagreements, assumptions, limitations, and unresolved questions.

## Quality expectations
- Cite the paper and supplementary sources.
- Explain how code and analyses can be reproduced.
- Add tests for non-trivial code.
- Do not commit passwords, API keys, private data, or copyrighted full-text papers.
- Distinguish the paper's original method from your own implementation or extension.
"""

CONCEPTUAL = """# Conceptual Questions (10)

**Owner:** Researcher A

For each question, provide an evidence-based answer and cite page, figure, table, or section where possible. Link relevant computational questions.

## Questions
For Q1–Q10, complete:
- **Question:**
- **Answer:**
- **Evidence from paper:**
- **Related computational question(s):** Q1–Q12
"""

COMPUTATIONAL = """# Computational / Coding Questions (12)

**Owner:** Researcher B

All 12 questions must be connected to one or more of Researcher A's conceptual questions. Focus on mathematical formulation, data structures, algorithms, implementation, complexity, testing, evaluation, and reproducibility.

## Questions
For Q1–Q12, complete:
- **Question:**
- **Linked conceptual question(s):** Q1–Q10
- **Computational reasoning / derivation:**
- **Implementation or pseudocode:**
- **Test / expected output:**
- **Interpretation and limitations:**
"""

DISCUSSION = """# Joint Discussion

- **Date:**
- **Contributors:**
- **Paper:**

## Discussion notes
- Key conceptual points:
- Key computational points:
- Assumptions challenged:
- Alternative methods considered:
- Points of agreement:
- Open disagreements:

## Final takeaways
- What did we learn?
- What can we reproduce?
- What remains uncertain?
- What should we investigate next?
"""

VALIDATOR = """#!/usr/bin/env python3
\"\"\"Validate the expected five-domain, 100-paper folder structure.\"\"\"
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
DOMAIN_ROOT = ROOT / "Domain"

def main():
    errors = []
    if not DOMAIN_ROOT.is_dir():
        errors.append("Missing Domain/ directory.")
    else:
        domains = sorted(p for p in DOMAIN_ROOT.iterdir() if p.is_dir())
        if len(domains) != 5:
            errors.append(f"Expected 5 domain folders; found {len(domains)}.")
        total = 0
        for domain in domains:
            if not (domain / "README.md").is_file():
                errors.append(f"Missing {domain.relative_to(ROOT)}/README.md")
            papers = sorted(p for p in domain.iterdir()
                            if p.is_dir() and p.name.startswith("paper-"))
            total += len(papers)
            if len(papers) != 20:
                errors.append(f"{domain.name}: expected 20 paper folders; found {len(papers)}.")
            for paper in papers:
                for filename in ("README.md", "conceptual-questions.md",
                                 "computational-questions.md", "discussion.md"):
                    if not (paper / filename).is_file():
                        errors.append(f"Missing {paper.relative_to(ROOT)}/{filename}")
        if total != 100:
            errors.append(f"Expected 100 paper folders total; found {total}.")
    if errors:
        print("Repository validation failed:")
        for error in errors:
            print("-", error)
        return 1
    print("Validation passed: 5 domains and 100 paper workspaces found.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
"""

CI = """name: Repository CI

on:
  push:
    branches: [main]
  pull_request:
    branches: [main]

jobs:
  validate:
    runs-on: ubuntu-latest
    steps:
      - name: Check out repository
        uses: actions/checkout@v4
      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: "3.12"
      - name: Validate domain and paper structure
        run: python scripts/validate_structure.py
      - name: Check Python syntax
        run: python -m compileall -q scripts
"""

def write_file(path: Path, content: str, overwrite: bool):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        print(f"Kept existing: {path}")
        return
    path.write_text(content.rstrip() + "\\n", encoding="utf-8")
    print(f"Wrote: {path}")

def paper_readme(paper_id: int, domain_name: str) -> str:
    return f"""# Paper {paper_id:03d}

## Paper information
- **Title:**
- **Authors:**
- **Year:**
- **Journal / conference:**
- **DOI / URL:**
- **Domain:** {domain_name}

## Research objective
- What computational problem does the paper address?
- What does it predict, explain, reconstruct, or optimise?

## Method overview
- **Input data:**
- **Network representation:**
- **Mathematical formulation:**
- **Algorithm / model:**
- **Baselines:**
- **Evaluation metrics:**
- **Main findings:**

## Completion checklist
- [ ] Both contributors read the paper
- [ ] 10 conceptual questions drafted
- [ ] 12 linked computational questions drafted
- [ ] Questions reviewed by both contributors
- [ ] Joint discussion completed
- [ ] Implementation completed (if applicable)
- [ ] Results and limitations documented
- [ ] Pull Request reviewed and merged

## Critical assessment
1. What are the method's key assumptions?
2. What are its computational limitations?
3. Is the evaluation sufficient to support the claims?
4. What could be reproduced or independently tested?
5. What is one defensible extension?
"""

def create_project(repo: Path, overwrite: bool = False):
    repo.mkdir(parents=True, exist_ok=True)
    write_file(repo / "README.md", ROOT_README, overwrite)
    write_file(repo / "CONTRIBUTING.md", CONTRIBUTING, overwrite)
    write_file(repo / "scripts/validate_structure.py", VALIDATOR, overwrite)
    write_file(repo / ".github/workflows/ci.yml", CI, overwrite)

    domain_root = repo / "Domain"
    table = "| Domain | Paper range | Count | Computational focus |\\n|---|---:|---:|---|\\n"
    table += "\\n".join(
        f"| {n}. {name} | {start}–{end} | 20 | {focus} |"
        for n, slug, name, start, end, focus in DOMAINS
    )
    domain_index = (
        "# Domains of Computational Network Biology\\n\\n"
        "This directory contains **5 domains and 100 papers**. Each domain has 20 paper workspaces.\\n\\n"
        + table +
        "\\n\\nEach paper workspace contains templates for 10 conceptual questions, "
        "12 linked computational/coding questions, joint discussion, implementation, and results.\\n"
    )
    write_file(domain_root / "README.md", domain_index, overwrite)

    for n, slug, name, start, end, focus in DOMAINS:
        ddir = domain_root / slug
        d_readme = (
            f"# Domain {n}: {name}\\n\\n"
            f"- **Paper range:** {start}–{end}\\n"
            f"- **Number of papers:** 20\\n"
            f"- **Computational focus:** {focus}\\n\\n"
            "Complete the 20 paper folders in numerical order where practical. "
            "Connect every computational question to the conceptual question(s) it develops.\\n"
        )
        write_file(ddir / "README.md", d_readme, overwrite)
        for paper_id in range(start, end + 1):
            pdir = ddir / f"paper-{paper_id:03d}"
            write_file(pdir / "README.md", paper_readme(paper_id, name), overwrite)
            write_file(pdir / "conceptual-questions.md", CONCEPTUAL, overwrite)
            write_file(pdir / "computational-questions.md", COMPUTATIONAL, overwrite)
            write_file(pdir / "discussion.md", DISCUSSION, overwrite)
            (pdir / "implementation").mkdir(parents=True, exist_ok=True)
            (pdir / "results").mkdir(parents=True, exist_ok=True)

    print("\\nDone. Review the generated files, then commit and push to GitHub.")

def main():
    parser = argparse.ArgumentParser(description="Create the Computational Network Biology project structure.")
    parser.add_argument("--repo", default=".", help="Cloned repository root (default: current directory).")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite generated documentation files.")
    args = parser.parse_args()
    create_project(Path(args.repo).resolve(), args.overwrite)

if __name__ == "__main__":
    main()
