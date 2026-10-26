# Endometriosis transcriptomic reanalysis

Source audit and sample-level reanalysis of published endometriosis transcriptomic features. Version 1.0.0; creator: [Huafeng Yu](https://orcid.org/0009-0005-9020-1272), Zhujiang Hospital, Southern Medical University.

This repository accompanies a manuscript prepared for BMC Medical Genomics. It reproduces fixed gene-association analyses, HC3 sensitivity analyses, and four figures from the supplied processed inputs. It does not refit published diagnostic classifiers or estimate AUC. The manuscript is not yet published.

## Reproduce

Clone this repository, or download and extract `endometriosis-transcriptomic-reanalysis-v1.0.0-full.zip` from the [v1.0.0 release](https://github.com/yu0598/endometriosis-transcriptomic-reanalysis/releases/tag/v1.0.0). Both include the processed matrices and numeric reference files. Install the exact R packages listed in [SOFTWARE.md](SOFTWARE.md), then run from the repository root:

```sh
python3 scripts/reproduce.py
```

Python uses only its standard library. Rscript must be on PATH; alternatively use `python3 scripts/reproduce.py --rscript /path/to/Rscript`. Set `ENDOMETRIAL_R_LIBRARY` if the matching packages are in a separate installed library. The runner creates a relative temporary directory inside this repository and stops if any step fails. No automatic downloads or package changes occur.

The run writes `results/` and `qa/`: nine fitted models; all 3,241 fixed target rows; primary and secondary HC3 tables (471 and 1,428 rows); and PDF, SVG, 600-dpi PNG/TIFF, and preview versions of four figures. Statistical results are compared against frozen references using unchanged absolute 2e-10 / relative 2e-9 tolerances. Figure scripts redraw the complete, frozen `figure_source/` tables; they do not select a new subset from the fitted results. See [docs/VALIDATION.md](docs/VALIDATION.md).

Individual steps, from the repository root:

```sh
mkdir -p runtime/tmp
TMPDIR=runtime/tmp Rscript --vanilla scripts/run_all.R
TMPDIR=runtime/tmp python3 scripts/compare_results.py
TMPDIR=runtime/tmp Rscript --vanilla scripts/run_hc3.R
python3 tests/check_hc3.py
TMPDIR=runtime/tmp Rscript --vanilla scripts/run_figures.R
```

## Study boundaries

- GSE51981: primary 144-sample, nine-cell model; fixed U = 157. Mean, median, 148-sample additive and identical global-QC sensitivity are retained. The secondary contrast family contains 471 tests.
- GSE120103: primary raw 34-sample model; fixed U2 = 119. All six prespecified processing/sample modes retain the full family, including unmeasured entries assigned P = 1. Mean, interaction and same-direction joint tests retain their original definitions.
- The primary same-direction joint test has 9 BH119 and 3 BY119 rejections. The three BY119 features are FOS, SYK and PDIA4. These are association results, not a validated diagnostic panel.
- The CEL audit concerns 61 records representing 42 distinct measurements and 19 exact cross-series pairs. Measurement identity is not patient identity or proof of leakage in a particular published classifier.
- Scanner/FE imbalance, the two incomplete raw FE files, and the unresolved relationship of the source-processed export to its baseline remain limitations.

## Contents and provenance

`data/`: exact processed matrices, frozen gene families, clinical and technical metadata. `reference/`: target results and background fitted coefficients used for validation. `tables/`: authoritative study CSV tables. `figure_source/`: complete figure tables. `scripts/`: portable entrypoints. `provenance/` and `optional_raw/`: historical methods, audits, upstream retrieval records and adaptation evidence.

Historical scripts are archival recipes requiring the larger original project and are not portable entrypoints. This public edition redacts personal absolute paths in historical text records; `provenance/original_processed_package_manifest.csv` records the earlier package hashes. No original statistical model, threshold or gene-family definition was changed. The current required-input hashes are in `SEMANTIC_INPUT_MANIFEST.csv`; the complete public-file inventory is `PUBLIC_FILE_MANIFEST.csv`.

See [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md), [DATA_LICENSE.md](DATA_LICENSE.md), [REFERENCES.md](REFERENCES.md), [optional_raw/README.md](optional_raw/README.md), and [CITATION.cff](CITATION.cff). MIT applies to original code; third-party data and identifiers retain their source terms. Full raw archives, journal full texts, internal review materials, manuscript drafts and installed runtimes are excluded.

Repository: https://github.com/yu0598/endometriosis-transcriptomic-reanalysis . Cite the fixed version and its archival DOI when available; no DOI is asserted until the archive is published.
