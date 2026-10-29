# Run every target from the repository root.
PY ?= python
RS ?= Rscript

.PHONY: figures tables genes panel infercnv ikarus copykat spikein preprint preprint-bib preprints supplement

figures:            ## Fig 1-4 from results/ (no raw data needed)
	$(PY) scripts/fig1.py
	$(PY) scripts/fig2.py
	$(PY) scripts/fig3.py
	$(PY) scripts/fig4.py

tables:             ## manuscript tables + derived numbers from results/
	$(PY) scripts/make_tables.py

genes:              ## Ensembl coordinates -> gene_positions.csv, loci_*.csv
	$(PY) scripts/gene_table.py

panel:              ## download + subsample CELLxGENE datasets -> panel/*.h5ad
	$(PY) scripts/assemble.py

infercnv: genes panel    ## inferCNV scores + decision rules -> calls/
	$(PY) scripts/run_infercnv.py

ikarus: panel       ## needs pyshim/ and ikarus_assets/ (see README)
	$(PY) scripts/run_ikarus.py

copykat: panel      ## export to Matrix Market, then CopyKAT in R
	$(PY) scripts/export_for_r.py
	$(RS) R/run_copykat_all.R

spikein: genes panel     ## needs chrom_arms.csv (see README)
	$(PY) scripts/build_spikein.py

supplement:         ## supplemental material zip for the preprint server
	$(PY) scripts/make_supplement.py

TEXBIN ?= /Users/jubayer/Library/TinyTeX/bin/universal-darwin

# pdflatex -> bibtex -> pdflatex x2 (refs, lastpage, forward refs to Tables S1-S6).
# Builds in $(1)/ next to the source; only the .tex, .bib and .pdf remain there.
define latex_build
	@export PATH="$(TEXBIN):$$PATH" && \
	cd $(1) && \
	pdflatex -interaction=nonstopmode $(2).tex >/dev/null && \
	bibtex $(2) >/dev/null && \
	pdflatex -interaction=nonstopmode $(2).tex >/dev/null && \
	pdflatex -interaction=nonstopmode $(2).tex >/dev/null && \
	if grep -qiE "undefined (reference|citation)" $(2).log; then \
	  echo "WARNING: unresolved references remain in $(2).pdf -- do not submit"; \
	else echo "built $(1)/$(2).pdf"; fi && \
	rm -f $(2).aux $(2).log $(2).out $(2).bbl $(2).blg
endef

preprint: tables      ## build the bioRxiv preprint PDF -> preprint/
	$(call latex_build,preprint,scCNV-FP-Evaluation)

preprint-bib: tables  ## journal variant -> journal/ (source kept locally, not in this repo)
	$(call latex_build,journal,scCNV-FP-Evaluation_BiB)

preprints: preprint preprint-bib  ## build both manuscript PDFs
