# External data inventory

Large source datasets are not copied into Git. Machine-specific locations belong
in `config/project_paths.local.yaml`, which is ignored. The versioned inventory is
`provenance/external_data_manifest.csv`.

Required source material:

- GSE102556 raw counts and GRCh38.p13 NCBI gene annotation;
- GEO sample metadata (a local clean export may be used, otherwise the historical R workflow used `GEOquery::getGEO("GSE102556")`);
- Allen Human Brain Atlas microarray donor data for the AHBA preprocessing branch;
- Hawrylycz et al. (2015) Supplementary Table 6 for published AHBA modules;
- MSigDB/GO/KEGG/Reactome annotation resources where a registered enrichment stage requires them.

Licensing and redistribution terms were not assumed. The helper script only
downloads sources with a confirmed public URL; restricted resources remain
manual inputs. Checksums, not copies, establish identity.
