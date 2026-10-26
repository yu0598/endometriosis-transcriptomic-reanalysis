# Software and installation

Observed reference environment: R4.6.1 (2026-06-24), limma3.68.5, statmod1.5.2, jsonlite2.0.0, digest0.6.39; macOS/aarch64 and BLAS/LAPACK details are in provenance/sessions/. Platform binaries are deliberately excluded. The fixed Entrez-to-symbol export comes from org.Hs.eg.db3.23.1 (EG source2026-Mar18); that package/database is not required to rerun the supplied matrices because the exact canonical membership table is included.

Install packages into your own platform's R library. In R4.6:

```r
install.packages(c("BiocManager", "statmod", "jsonlite", "digest"))
BiocManager::install(version="3.23", ask=FALSE, update=FALSE)
BiocManager::install("limma", ask=FALSE, update=FALSE)
sapply(c("limma","statmod","jsonlite","digest"), packageVersion)
```

The portable entrypoints require the exact four reference package versions and stop on mismatches. Repository versions can change: obtain the matching source release from the official package archive if the current repository no longer supplies it; compiling from source may require your platform's R build tools. No automatic downgrade or package installation occurs during a statistical run. The install commands are documented from official repositories; this task tested the rerun using the already installed matching library, not a fresh internet installation.

Official sources checked2026-10-26: [limma/Bioconductor](https://bioconductor.org/packages/limma/), [statmod/CRAN](https://cran.r-project.org/package=statmod), [jsonlite/CRAN](https://cran.r-project.org/package=jsonlite), [digest/CRAN](https://cran.r-project.org/package=digest). The independent historical Python checks used NumPy2.5.4, SciPy1.18.1 and Statsmodels0.15.0; those packages are optional here because only the R statistical rerun and standard-library comparator are portable entrypoints.

## Figure dependencies

The figure scripts additionally use ggplot2 4.0.3, cowplot 1.2.0, ragg 1.5.2 and svglite 2.2.2. These versions were tested locally; the plotting scripts do not enforce exact versions. Fonts may differ between platforms. Cairo PDF is preferred when supported by the R build; SVG is always provided as editable text. Statistical entrypoints enforce the four reference package versions above.
