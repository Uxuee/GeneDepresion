# External data acquisition and placement

Large, licensed, or source-controlled datasets are not redistributed through
this Git repository. Their identities and expected destinations are recorded in
`provenance/external_data_manifest.csv`.

## Configure local paths

Copy `config/project_paths.example.yaml` to
`config/project_paths.local.yaml`, then set:

- `raw_data_root` for GSE102556 counts, annotation, and optional metadata;
- `external_resource_root` for AHBA and published supplementary resources;
- `output_root` for generated canonical outputs;
- `r_library_path` if using a separate R library.

The local configuration file is ignored by Git.

## Required inputs

| Resource | Source/accession | Expected local destination | Identity/acquisition |
|---|---|---|---|
| GSE102556 counts | NCBI GEO `GSE102556` | `${raw_data_root}/GSE102556_raw_counts_GRCh38.p13_NCBI.tsv` | Manual GEO acquisition; expected SHA-256 is in the manifest |
| GSE102556 metadata | NCBI GEO `GSE102556` | `${raw_data_root}/GSE102556_metadata_clean.csv` | `GEOquery::getGEO("GSE102556")` or a documented clean export |
| GRCh38.p13 annotation | NCBI | `${raw_data_root}/Human.GRCh38.p13.annot.tsv` | Manual acquisition; expected SHA-256 is in the manifest |
| AHBA donor microarray files | Allen Human Brain Atlas | `${external_resource_root}/AHBA_microarray/` | Obtain under Allen Institute terms; exact historical donor-file manifest remains incomplete |
| Hawrylycz 2015 Table 6 | Published supplement | `${external_resource_root}/Hawrylycz_2015_Supplementary_Table_6.csv` | Public helper available; expected SHA-256 is in the manifest |
| GO/Reactome/MSigDB/KEGG resources | Package/resource providers | R library/resource cache | Restore only the versions evidenced by the relevant stage environment |

Download the URL-confirmed published module table without overwriting an
existing copy:

```text
python code/download_public_external_data.py hawrylycz_table6 --destination /path/to/external/resources
```

Check local availability and known checksums:

```text
python code/run_canonical_pipeline.py --check-external-data
```

Missing required data produce an actionable list of destinations and setup
instructions. Package-managed resources are reported separately and are not
mistaken for missing files.

## Licensing and redistribution

The manifest records the known licensing/redistribution status. `UNKNOWN` means
that redistribution permission was not established—it does not mean that the
file is public-domain. Do not commit raw data merely to make a clone
self-contained. Acquire each resource from its provider and verify checksums
where recorded.
