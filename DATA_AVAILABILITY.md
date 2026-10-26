# Data and code availability

Code and the processed-input reproducibility package are available at https://github.com/yu0598/endometriosis-transcriptomic-reanalysis . The fixed v1.0.0 software archive is identified by https://doi.org/10.5281/zenodo.23670029 .

Included: exact processed gene matrices, fixed U157/U2=119, sample masks and metadata, frozen canonical Entrez membership, target/coefficient references, study CSVs, figure source tables, historical methods/audits and hash manifests. The raw-arm RDS files contain numeric matrices, not platform executables.

Underlying public deposits: [GSE51981](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51981), [GSE120103](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120103), [GSE7305](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7305), [GSE7307](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7307), [GPL570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570), [GPL6480](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6480). File URLs, sizes and original hashes are in optional_raw/source_manifests/. Accessions identify the original source deposits, not a new deposit made by this study.

Not included: full raw archives/FE gzip files, installed runtimes, platform libraries, the complete original annotation database, copyrighted paper collections, manuscript drafts and internal author reviews. No patient crosswalk or unreported clinical covariate was inferred. The tested reproducibility claim begins at the included processed inputs after installation of matching dependencies. Raw normalization and upstream raw reconstruction were not rerun for this release. The source-processed GSE120103 baseline/export relationship remains unresolved.
