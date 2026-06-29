# AHBA Data Directory

Raw Allen Human Brain Atlas microarray downloads should be stored here locally
but are intentionally not committed to git.

Use:

```bash
python code/ahba_modules/00_fetch_ahba_data.py
```

The script first checks this folder for existing AHBA files. If they are missing
and `abagen` is installed with internet access, it will try:

```python
abagen.fetch_microarray(donors="all", data_dir="data/ahba")
```

If automatic download fails, download AHBA microarray files manually into this
folder, keeping donor-level files such as `MicroarrayExpression.csv`,
`Probes.csv`, and `SampleAnnot.csv` in donor subfolders when possible.
