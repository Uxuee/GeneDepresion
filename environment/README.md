# Environment provenance

The repository contains different levels of environment recovery. Do not treat
the older root `environment.yml` as a fully pinned publication environment.

## Python

Canonical Python validation uses the standard-library `unittest` runner plus
the numerical packages imported by the project scripts (`pandas`, `numpy`,
`scipy`, and scikit-learn where required). The exact original end-to-end Python
environment is not fully recovered.

On the audited Windows machine, the bundled runtime was blocked by a local ACL,
so validation used an ACL-safe temporary Python copy. That temporary path is
machine-local evidence, not a dependency to commit or assume on a new machine.
A normal supported Python installation may be used instead.

Run the test suite with:

```text
python -m unittest discover -s tests -p "test*.py" -v
```

## Phase‑3B R/WGCNA environment

Phase 3B captured a genuine reconstructed-current R environment under
`environment/phase3b/`:

- `renv.lock` — R 4.6.1, Bioconductor 3.23, and the exact recorded package set;
- `package_versions.csv` — compact package-version inventory;
- `sessionInfo.txt` — the session used for the approved Phase‑3B validation;
- `install_phase3b_packages.R` and `snapshot_environment.R` — reconstruction
  helpers.

The recorded session includes WGCNA 1.74, dynamicTreeCut 1.63-1,
clusterProfiler 4.20.0, ReactomePA 1.56.0, org.Hs.eg.db 3.23.1, GO.db 3.23.1,
and msigdbr 26.1.1. These are reconstructed-current versions used by Phase 3B;
they are not claims about unknown historical runs.

Typical restoration from the repository root is:

```text
Rscript -e "renv::restore(lockfile='environment/phase3b/renv.lock', prompt=FALSE)"
Rscript audit/validate_phase3b_wgcna.R
```

The committed `sessionInfo.txt` is canonical evidence. A local
`sessionInfo_after_install.txt`, if produced during troubleshooting, is not a
canonical artifact and should remain untracked.

## Historical/partly recovered R stages

Older gene-model and enrichment stages retain recorded versions where known,
but their complete historical dependency graphs were not reconstructed. The
canonical registry marks stages with insufficient provenance as historical,
pending recovery, or excluded rather than inventing package versions.
