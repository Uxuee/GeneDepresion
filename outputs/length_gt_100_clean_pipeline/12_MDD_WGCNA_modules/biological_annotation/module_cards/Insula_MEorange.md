# Insula MEorange

- Module size: 53
- Diagnosis beta/p/FDR: -0.0939 / 0.0324 / 0.469
- Tier2/ML enrichment summary: same_region_Tier2_genes overlap=2 FDR=0.219; same_region_high_plus_moderate_ML_genes overlap=1 FDR=0.944; all_Tier2_genes overlap=2 FDR=0.961; high_plus_moderate_ML_genes overlap=1 FDR=1; high_confidence_ML_genes overlap=0 FDR=1
- AHBA overlap summary: best M2 FDR=1.08e-04; M2 overlap=6 FDR=1.08e-04; M26 overlap=2 FDR=1; M28 overlap=0 FDR=1; M6 overlap=4 FDR=0.363; M7 overlap=1 FDR=1; Top Hawrylycz enrichment terms point to hsa-miR-448:PITA_TOP; central nervous system neuron development; hsa-miR-548l:PITA_TOP; hsa-miR-606:PITA_TOP; hsa-miR-548n:PITA_TOP. Categories represented: MicroRNA, GO: Biological Process.
- Top hub genes: OLFML2B; HS3ST4; GRIK3; GALNT14; PTPRU; SMYD1; CCN4; SHB; NR4A2; TMEM233
- Top GO BP terms: Commitment of neuronal cell to specific neuron type in forebrain; Forebrain neuron fate commitment; Glial cell derived neurotrophic factor receptor signaling pathway
- Top GO CC terms: External encapsulating structure; Node of ranvier; Axon
- Top pathway terms: Apical surface; Sig insulin receptor pathway in cardiac myocytes; Sig pip3 signaling in cardiac myoctes
- Cell-type marker summary: no clear marker signal

## Careful Interpretation

Insula MEorange is a prioritized exploratory module. Its annotation is consistent with neuronal/synaptic biology based on enrichment, AHBA context, and hub-gene patterns. However, its module eigengene did not survive BH-FDR correction for Diagnosis, so it should be treated as a follow-up candidate rather than a validated MDD-associated module.

## Caution

Exploratory module annotation only. No module eigengene survived Diagnosis BH-FDR < 0.10 in the full Step6 run.
