# Environment provenance

The repository does not yet have one fully pinned, end-to-end environment lock.
No package versions were invented during Phase 3A.

## Recovered exact information

- Historical WGCNA annotation manifest: R 4.4.2 and `msigdbr` 26.1.0; `org.Hs.eg.db`, `clusterProfiler`, and `ReactomePA` were recorded as unavailable for that run.
- Historical dream run documentation: R 4.4.2, `variancePartition` 1.36.3, `lmerTest` 3.2.1, and `lme4` 1.1.37.
- Python scoring tests are runnable with the temporary recovered runtime at `C:/Users/palom/AppData/Local/Temp/codex-python-ml-audit/python.exe`. That runtime is not a project artifact and must not be committed.

`environment.yml` is a historical, incomplete environment description rather
than a publication-grade lock. A future controlled rerun should capture R
`sessionInfo()`, Python package versions, operating system, and database release
identifiers in a new run manifest. Until then, stages with insufficient
environment recovery are marked `pending_recovery` or excluded.

## Stage requirements

| Area | Runtime/dependencies evidenced in repository | Version status |
|---|---|---|
| Python scoring/downstream | Python with pandas, numpy, scipy, scikit-learn, matplotlib/seaborn and network utilities used by the scripts | exact full historical set unknown; temporary test runtime currently works |
| Python tests | standard-library `unittest`; scoring tests import the project scorer and installed numerical stack | currently tested with temporary runtime; not locked |
| R gene models | R/Bioconductor, limma/edgeR, variancePartition/dream where applicable | selected exact historical versions above; full dependency graph unknown |
| WGCNA | R, WGCNA and dynamic tree-cut dependencies | saved run manifest is the environment reference; complete lock unavailable |
| Enrichment | R/Python enrichment scripts plus MSigDB/GO/KEGG/Reactome resources | resource snapshots partly unknown |
| AHBA | Python data processing/ML and R WGCNA for de novo branch | full historical package set unknown |

`Rscript` was not available in the audited command environment, so no new
`sessionInfo()` was captured. A requirements or `renv.lock` file was not created
because doing so without a fully usable environment would falsely imply exact
reproducibility.
