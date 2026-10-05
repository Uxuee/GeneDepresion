# ML–RRA experimental environment

Python 3.12.14 was used with the packages pinned in `requirements.txt`.

R 4.6.1 was used with `RobustRankAggreg` 1.2.1. The installed R library and copied Python runtime are machine-local and ignored. Restore them with:

```powershell
python -m pip install -r environment/ml_rra/requirements.txt
Rscript -e "dir.create('environment/ml_rra/library', recursive=TRUE, showWarnings=FALSE); install.packages('RobustRankAggreg', repos='https://cloud.r-project.org', lib='environment/ml_rra/library')"
```

Set `R_LIBS_USER` to the absolute `environment/ml_rra/library` directory when running the RRA stage.
