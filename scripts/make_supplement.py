"""Build the supplemental-material zip for the preprint server.

Collects the result tables every manuscript number is read from, the analysis and
figure scripts, the figure files and the LaTeX source. Run from the repository root:

    python scripts/make_supplement.py
"""
import os
import zipfile

OUT = "scCNV-FP-Evaluation_supplement.zip"

PATTERNS = [
    ("results", (".csv", ".json", ".md")),
    ("scripts", (".py",)),
    ("R", (".R",)),
    ("figures", (".pdf", ".png")),
    ("preprint", (".tex", ".bib")),
    ("tables", (".tex",)),
    ("config", (".json",)),
    ("metadata", (".json",)),
]
# docs/ is deliberately excluded: it holds internal drafting and review records
# (docs/peer_review_report.md, docs/scientific_review_*.md, submission checklists)
# that are not reader-facing supplementary material.
SINGLES = ["README.md", "LICENSE", "Makefile", "environment.yml",
           "requirements-figures.txt", "calls/TS_Lung_infercnv.parquet"]


def main():
    if not os.path.isdir("results"):
        raise SystemExit("run from the repository root")
    n = 0
    with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as z:
        for d, exts in PATTERNS:
            if not os.path.isdir(d):
                continue
            for f in sorted(os.listdir(d)):
                p = os.path.join(d, f)
                if "_BiB" in f:
                    continue  # journal-variant source, not part of the preprint supplement
                if os.path.isfile(p) and f.endswith(exts) and not f.startswith("."):
                    z.write(p, p)
                    n += 1
        for p in SINGLES:
            if os.path.isfile(p):
                z.write(p, p)
                n += 1
    print(f"wrote {OUT} with {n} entries "
          f"({os.path.getsize(OUT) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
