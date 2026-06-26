# Step 9: ML v1 Tier2_fixed Only

This clean branch filters existing ML v1 outputs to feature_set == Tier2_fixed.
Tier4_top25/50/100/250/500 and strict_univariate_top100 are excluded.
ML v1 uses Tier 2 gene expression features to separate MDD vs Control samples by cross-validation within the same dataset.
Because Tier 2 genes were selected using the full dataset, performance metrics are exploratory. The main purpose is gene prioritization.
