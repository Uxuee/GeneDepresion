# Step 9: Extended ML v1 Prioritization

This clean branch reuses existing ML v1 outputs and keeps Tier2_fixed plus Tier4_top-k feature sets.
Tier 2 genes define the main exploratory candidate pool. Tier 4 ranked signatures are used only as an additional ML-support layer to test whether candidate-gene prioritization is stable across feature-set definitions.
strict_univariate_top100 is retained only as a reference/sanity branch, not as the main candidate source.
ML v1 is interpreted as gene prioritization, not validated diagnostic prediction. Classification performance is descriptive because feature sets were derived from the same dataset.
