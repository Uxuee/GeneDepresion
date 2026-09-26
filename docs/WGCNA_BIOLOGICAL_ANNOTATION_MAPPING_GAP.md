# WGCNA biological-annotation mapping gap

The exact historical annotator is `code/Step7_Annotate_MDD_WGCNA_Modules.R`.
Its manifest records R 4.4.2 and `msigdbr` 26.1.0, with `org.Hs.eg.db`,
`clusterProfiler`, and `ReactomePA` unavailable.

The saved mapping report contains 4,267 rows: 1,796 mapped (42.1%) and 2,471
unmapped (57.9%). Approximate historical coverage varied strongly by analysis
unit: pooled Frontal Cortex 25.2%, Hippocampal Formation 44.0%, Insula 64.8%,
and Ventral Striatum 70.4%. The script checks package availability, but the
saved mapping report is based on the `msigdbr` mapping join; it does not provide
evidence that an `org.Hs.eg.db` fallback actually mapped these rows.

`org.Hs.eg.db` is not currently testable in a reproducible R environment because
`Rscript` is unavailable in the audited shell and no complete R lock exists.
Accordingly, no reannotation was run. A future correction must go into a new
versioned directory, record input/script/output hashes and `sessionInfo()`, and
compare coverage, newly mapped genes, pathway/theme changes, and manuscript
interpretation. It must not change WGCNA modules, eigengenes, enrichment, or
structural-core definitions.

Until then, historical biological themes are limited annotation/context rather
than primary statistical evidence. Any manuscript claim that relies directly on
the incomplete mapping must carry this limitation or be deferred.
