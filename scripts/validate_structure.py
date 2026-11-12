#!/usr/bin/env python3
"""Validate the expected five-domain, 100-paper folder structure."""
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
    sys.exit(main())\n