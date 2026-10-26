# Release validation

On 2026-10-26, all nine statistical model fits were rerun from the shipped processed inputs using R 4.6.1 and the pinned dependencies in SOFTWARE.md. The comparison checked 3,241 target rows and the full fitted-coefficient backgrounds: 997,535 numeric cells, maximum absolute error 0. All row identities, metadata and missingness checks passed.

The separate portable HC3 runner recomputed 471 primary and 1,428 secondary rows. The standard-library test compared 4,695 and 14,060 numeric fields respectively against the frozen independent references. Maximum absolute differences were 1.0671463712697005e-12 and 8.859579736508749e-14, below the unchanged absolute 2e-10 / relative 2e-9 tolerance. This validates the fixed HC3 output, not a new primary analysis.

All four R plotting entrypoints ran successfully; all four previews were visually inspected for readable labels and unclipped content. Figures redraw complete frozen source tables and retain the original geometry, estimands, selection annotations and counts. SVG text remains editable. PDF uses Cairo when available and base R PDF otherwise; raster outputs use ragg. Source-adaptation details are in metadata/figure_adaptations.json.

Machine-readable reports are in docs/validation/. The earlier independent audits remain in provenance/audits/. Their one-ULP source-parser and near-flat background warnings are retained. Historical reports describe their original scope; they are not proof that every historical script can be rerun with this compact package.

Limits: one local macOS/aarch64 environment was exercised; a fresh internet dependency installation and cross-platform equality were not tested. Full raw reconstruction, prospective biological validation and classifier performance are outside this release.

A separate read-only code review re-solved OLS for all 1,852 estimable HC3 targets using the full sandwich covariance matrix; maximum SE difference was 1.166e-15. Mutation checks rejected perturbed SE, missing-test P and adjusted P. The test now also asserts exact P = 1 and empty effect/uncertainty fields for every non-estimable target.

Static figure checks passed without FAIL findings. Remaining warnings are heuristic: the validator does not invoke the installed R parser, the workflow grid has no statistical uncertainty panel, and widths are stored as variables rather than device literals. The actual R scripts executed successfully, all figures use 183 mm final width, and numeric panels retain their stated model intervals.
