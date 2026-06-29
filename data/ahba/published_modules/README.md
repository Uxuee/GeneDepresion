# Published AHBA Canonical Modules

Place the published Hawrylycz et al. 2015 AHBA canonical module assignment table here before running the published-module enrichment branch.

Preferred file: Supplementary Table 6 from Hawrylycz et al. 2015, *Canonical genetic signatures of the adult human brain*. This table should contain gene symbols and the published adult human brain consensus module/signature assignments. Eigengene/module-membership correlations are used when available.

Supported formats:

- `.csv`
- `.tsv` or `.txt`
- `.xlsx` or `.xls`

Suggested filename examples:

```text
Hawrylycz_2015_Supplementary_Table_6.xlsx
Hawrylycz_2015_Supplementary_Table_6.csv
```

Then run from the project root:

```bash
python code/ahba_modules/04_fetch_or_load_published_ahba_modules.py
python code/ahba_modules/05_test_published_ahba_module_enrichment.py
```

Raw downloaded AHBA data and supplementary module tables are not committed to git. This README is the only file in this folder intended for version control.


Official Supplementary Table 6 CSV URL used by the loader:

```text
https://static-content.springer.com/esm/art%3A10.1038%2Fnn.4171/MediaObjects/41593_2015_BFnn4171_MOESM97_ESM.csv
```
