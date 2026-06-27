# Step 10: Extended ML-Prioritized Tier 2 Genes

This step reuses the previous stable ML-prioritized Tier 2 output package.
Tier 2 genes define the main exploratory candidate pool. Tier 4 ranked signatures are used only as an additional ML-support layer to test whether candidate-gene prioritization is stable across feature-set definitions.
Tier 4 top-k feature sets are selected by absolute signed standardized log-p score, while the exported tables keep signed_score, abs_score, and direction so up-in-MDD and down-in-MDD genes remain distinguishable.
The final candidates are ML-prioritized Tier 2 genes, not Tier 4 genes.
The final prioritized genes should be interpreted as internally supported candidate genes, not validated biomarkers.
