# ML priority score code audit

> **Historical pre-correction audit.** Its 180-event v2 denominator and implementation findings describe the legacy scorer. Current method: `ML_PRIORITY_SCORE_FINAL.md`.

## Scope and bottom line

This is a read-only implementation audit of the gene-level ML-prioritisation score. It reconstructs what the current code actually computes; it does not rerun the ML pipeline or change any results.

**Overall verdict: the arithmetic implementation of the exported score and confidence thresholds is reproducible, but the meaning documented elsewhere is only partly aligned with the code.** The most consequential findings are:

1. The final v1 score uses only the fixed Tier2 feature set, although the existing method text says that Tier4 top-k and strict-univariate feature sets contribute.
2. The empirical denominators are **125 v1 events** and **180 v2 events per region**, not 450 v2 opportunities.
3. The v2 frequency rewards repeated inclusion across nested top-k lists; it is not normalised only over runs in which a gene was eligible.
4. Random-forest importances are treated as positive values in the sign-consistency calculation, so that field is not purely agreement among signed linear coefficients.
5. The majority baseline contributes to v2 selection frequency, despite not contributing to the five-family model-support count.
6. Cross-validated classification performance does not enter the gene score.
7. The final 5% statistical-evidence component is region-wise min-max scaling of saved gene-level mixed-model evidence; it is not called “statistical_model_evidence_scaled” in code.

These issues do not imply that the exported numbers were calculated incorrectly. They mean that some components measure something different from their labels or prose descriptions and should be described more precisely.

## Files inspected

Primary current scoring implementation:

- D:/genes/GeneDepresion/code/extract_stable_ML_prioritized_Tier2_genes.py

Historical scripts preserved in the pipeline snapshot:

- D:/genes/GeneDepresion/outputs/length_gt_100_clean_pipeline/code_snapshot_2026-06-27_tier4_score_export/code/run_ML_pilot_Tier2_length_gt_100.py
- D:/genes/GeneDepresion/outputs/length_gt_100_clean_pipeline/code_snapshot_2026-06-27_tier4_score_export/code/run_ML_strict_v2_length_gt_100.py

Saved descriptions and outputs:

- D:/genes/GeneDepresion/outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_extended/run_manifest.json
- D:/genes/GeneDepresion/outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_extended/final_gene_ranking_method.md
- D:/genes/GeneDepresion/outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_extended/ML_prioritized_Tier2_extended_all_ranked.csv

The raw fold-level input paths named by the manifest are no longer present in the current working tree. Therefore, the final score can be reconstructed exactly from the saved ranked table, and event rules can be reconstructed from code, but a gene-by-gene list of every original fold/model record cannot be independently re-read from retained raw files.

## End-to-end implementation map

~~~text
normalised expression and metadata
  |
  +-- v1 feature sets and repeated CV
  |     original run generated fold-level importance rows
  |     current final scoring filters those rows to Tier2_fixed only
  |
  +-- strict v2 nested feature selection
  |     feature selection and tuning occur inside each outer training split
  |     selected genes are written for every model/top-k/method/fold event
  |
  +-- annotate each saved row as selected or not selected
  |
  +-- aggregate selected events by region and Tier2 gene
  |
  +-- combine v1 and v2 summaries with gene-level statistical evidence
  |
  +-- calculate weighted score and assign high/moderate/low confidence
~~~

## 1. What one “event” means

The scoring script defines an event by the available combination of:

- analysis level;
- region;
- feature set;
- model;
- repeat;
- fold;
- outer repeat; and
- outer fold.

Only columns that exist in a particular input table are used. Within each event, saved features are ranked by absolute importance. The event-level row count is the number of feature rows actually written for that event, not necessarily the full gene universe.

The aggregation denominator is the number of distinct event combinations present anywhere in that region’s event table. It is **region-global**, not gene-specific. A gene absent from an event gets no selected count for that event, while the event remains in the denominator.

## 2. Models and feature sets

The five model families counted for model support are:

| Internal model name | Interpretation | Signed coefficient? |
|---|---|---:|
| logistic_l1 | L1-penalised logistic regression | yes |
| elastic_net_logistic | elastic-net logistic regression | yes |
| logistic_l2 | L2-penalised logistic regression | yes |
| linear_svm | linear support-vector classifier | yes |
| random_forest_exploratory | random forest | no; non-negative feature importance |

The historical v1 script constructed:

- Tier2_fixed;
- Tier4_top25, top50, top100, top250, and top500;
- strict_univariate_top100, selected within each training fold.

However, the current scoring entry point explicitly filters v1 fold-level rows to:

~~~text
feature_set == "Tier2_fixed"
~~~

before annotation and aggregation. Consequently, **the current exported v1 contribution is fixed-Tier2-only**. Tier4 top-k and v1 strict-univariate events do not enter the final score produced by this script.

The strict v2 run uses two train-only selection methods:

- ANOVA/F-test; and
- absolute standardised effect size.

Each method is evaluated at top-k values 25, 50, 100, 250, and 500.

## 3. How “selected” is defined

Selection is not identical across model types.

| Branch | Code rule used in final scoring |
|---|---|
| L1 logistic | absolute coefficient greater than 1e-10 |
| Elastic-net logistic | absolute coefficient greater than 1e-10 |
| L2 logistic | within the top 20% of saved event rows by absolute importance |
| Linear SVM | within the top 20% of saved event rows by absolute importance |
| Random forest | within the top 20% of saved event rows by absolute importance |
| Other model, including majority baseline in v2 | saved selected flag, then nonzero flag, then absolute importance greater than 1e-10, in that order of availability |

For the three ranked models, top 20% means:

~~~text
importance rank <= ceiling(number of saved rows in that event * 0.20)
~~~

Thus it is not correct to say that all non-zero L2 or SVM coefficients are counted as selected in the final aggregation. Their raw writers may save nearly all non-zero coefficients, but the final scoring layer selects only the top 20% of rows within an event.

Sparse-model support is a separate Boolean indicator: at least one non-zero L1 or elastic-net coefficient was observed.

### Majority baseline caveat

The strict-v2 writer records every train-selected feature for every model, including the majority baseline, with selected=true and no coefficient. The final annotation treats this as an “other model” selection. Therefore the baseline can increase the v2 selected-event count and v2 selection frequency. It is not among the five model families used for model-support count.

This asymmetry is an implementation concern: a no-feature classifier contributes to gene stability merely because the feature selector supplied that gene to the baseline event.

## 4. Selection-frequency denominators

The code calculates, within each region and source:

~~~text
selection frequency =
number of selected gene-event rows /
number of distinct events present for the region
~~~

It does not divide by the number of events in which that particular gene was inside a top-k input set.

### v1 denominator

The historical v1 fixed-Tier2 design used:

- five random seeds/repeats;
- five folds per seed; and
- five scored model families.

This gives:

~~~text
5 repeats * 5 folds * 5 model families = 125 events per region
~~~

The saved final table confirms this exactly: every v1 frequency multiplied by 125 is an integer.

### v2 denominator

With the executed settings represented in the final table, strict v2 used:

- two feature-selection methods;
- five top-k sizes;
- six models when the majority baseline is included; and
- three outer folds.

This gives:

~~~text
2 methods * 5 top-k values * 6 models * 3 outer folds = 180 events per region
~~~

The saved final table confirms this exactly: every v2 frequency multiplied by 180 is an integer, within floating-point precision.

A proposed denominator of 450 does not reproduce the saved frequencies. The script supports outer repeats, but the executed saved result is consistent with one outer repeat and three outer folds.

### Meaning of nested top-k frequency

A gene selected in the top 25 is normally also eligible in top 50, 100, 250, and 500. It can therefore accumulate repeated support across nested lists. Conversely, a gene absent from top 25 still has the top-25 events in its denominator. The v2 frequency is best described as **support across the complete pre-specified method-by-k-by-model-by-fold grid**, not conditional model stability after eligibility.

## 5. Importance, rank, and model support

Within each event:

1. absolute importance is calculated;
2. rows are ranked in descending order using minimum ranks for ties; and
3. selection is assigned using the model-specific rules above.

The aggregate mean and median absolute importance and mean rank are calculated over **selected rows only**, not all attempted events.

Raw coefficient sizes and random-forest importances are on different scales. Although their selected-row raw means are retained, the final weighted score does not directly add those heterogeneous magnitudes. It uses a regional transformation of mean importance rank:

~~~text
importance_rank_score =
clip(1 - (gene mean rank - 1) / maximum regional mean rank, 0, 1)
~~~

Because the denominator is the maximum rank rather than maximum rank minus one, the worst observed rank receives 1/max_rank rather than exactly zero.

When both sources are present, the code first averages the v1 and v2 source-level mean ranks with equal source weight, ignoring missing values. It does not weight them by their numbers of events.

Model-support count is the number of distinct families among the five named model families that selected the gene at least once across v1 or v2. One selected occurrence is sufficient. The score scales it by the largest observed model-support count in that region:

~~~text
model_support_scaled =
gene model-family count / regional maximum model-family count
~~~

This is not the fraction of all fits supported by each family; frequency captures the repeated-event component.

## 6. Sign consistency

For each source, the implementation:

1. keeps selected-row importance values that are non-missing and have absolute value greater than 1e-10;
2. counts positive and negative values; and
3. returns the larger of positive fraction and negative fraction.

Thus:

~~~text
sign consistency = max(P(value > 0), P(value < 0))
~~~

The score is 1 when all usable values share a sign and approaches 0.5 when signs split evenly. With both sources present, v1 and v2 source summaries are averaged equally.

### Important semantic limitation

The implementation does not restrict this calculation to signed linear models. Random-forest importance is always non-negative and is included as positive evidence. This can inflate positive-direction agreement and can influence the reported direction. Majority-baseline values are missing and are discarded.

Accordingly, the existing description “agreement across signed linear coefficients” is not an exact description of the implemented field. A clean signed-direction metric should exclude random forest and any other unsigned model before counting signs.

Missing sign consistency is filled with 0.5 in the weighted score. For the high-confidence threshold only, missing consistency is filled with 1 before testing the 0.70 requirement. That second behaviour is permissive and should be stated or changed in a future revision.

## 7. Combining v1 and v2

If strict-v2 output exists for a gene, combined selection frequency is:

~~~text
combined_selection_frequency =
0.60 * v2_selection_frequency +
0.40 * v1_selection_frequency
~~~

If no v2 result exists, the code uses the v1 frequency alone.

The v2-support component is binary:

- 1.0 if v2 selection frequency is greater than zero;
- 0.0 if v2 data exist but the gene has no v2 support; and
- 0.5 when no v2 source is available at all.

It does not encode the magnitude of v2 stability beyond the separate combined-frequency term.

## 8. Gene-level statistical evidence

The final score merges a saved gene-level field named lmm_evidence. Within each region it is min-max scaled:

~~~text
lmm_evidence_scaled =
(gene evidence - regional minimum) /
(regional maximum - regional minimum)
~~~

If a region has no range, the helper returns 0.5 for non-missing values. Missing values are also filled with 0.5 when the final score is built.

The code does not expose a field named statistical_model_evidence_scaled. “Statistical model evidence” is a reasonable plain-language description, but reports should name the underlying input and define precisely how lmm_evidence was originally constructed before interpreting this component.

## 9. Exact final score

The current scoring implementation is:

~~~text
ML priority score =
  0.40 * combined selection frequency
+ 0.15 * model-support scaled
+ 0.10 * sparse-model support
+ 0.15 * importance-rank score
+ 0.10 * sign-consistency score
+ 0.05 * regional min-max lmm evidence
+ 0.05 * v2-support score
~~~

The result is clipped to the interval 0 to 1.

Sparse-model support is 1 when either L1 or elastic net selected the gene at least once and 0 otherwise.

### Confidence labels

High confidence requires all of:

- score at least 0.70;
- combined selection frequency at least 0.50;
- at least two model families;
- sparse-model support; and
- sign consistency at least 0.70, with missing sign consistency treated as 1 for this test.

Moderate confidence requires all of:

- score at least 0.45;
- combined selection frequency at least 0.25; and
- at least one model family.

All other genes are low confidence. High is evaluated first.

There is no additional explicit sign-disagreement exclusion beyond the high-confidence sign threshold. Moderate confidence can therefore include inconsistent signs.

## 10. Worked real-gene reconstruction

The following row exists in the saved final ranking:

| Field | Value |
|---|---:|
| Region | Cingulate gyrus 25 / CG25 |
| Gene | LOC107985795 |
| Gene ID | 107985795 |
| v1 selection frequency | 0.864 |
| v2 selection frequency | 0.4777777778 |
| Combined selection frequency | 0.6322666667 |
| Model-family count | 5 |
| Sparse-model support | true |
| Combined mean importance rank | 8.3711240310 |
| Sign consistency | 1.0 |
| lmm_evidence | 3.8651490534 |
| v2 support | true |
| Exported score | 0.8311053330 |
| Exported label | high |

The frequencies imply:

~~~text
v1 selected count = 0.864 * 125 = 108 events
v2 selected count = 0.4777777778 * 180 = 86 events
~~~

Because the original raw fold tables are absent, the retained data cannot identify exactly which individual models, folds, and v2 top-k lists supplied those 108 and 86 selections. The auditable branch-level trace is:

| Source | Event grid used by final score | Denominator | Selected events |
|---|---|---:|---:|
| v1 | Tier2_fixed, five model families, five repeats, five folds | 125 | 108 |
| v2 | two selectors, five top-k values, six models including baseline, three outer folds | 180 | 86 |

### Component reconstruction

For CG25, the retained ranked table gives:

- regional maximum model count = 5;
- regional maximum combined mean rank = 33.25;
- regional minimum lmm_evidence = 1.6793093169; and
- regional maximum lmm_evidence = 6.1079366624.

Therefore:

~~~text
combined frequency
= 0.60(0.4777777778) + 0.40(0.864)
= 0.6322666667

model-support scaled
= 5 / 5
= 1

sparse support
= 1

importance-rank score
= 1 - (8.3711240310 - 1) / 33.25
= 0.7783120592

sign-consistency score
= 1

lmm-evidence scaled
= (3.8651490534 - 1.6793093169) /
   (6.1079366624 - 1.6793093169)
= 0.4935704827

v2-support score
= 1
~~~

Substitution into the score:

~~~text
0.40(0.6322666667)
+ 0.15(1)
+ 0.10(1)
+ 0.15(0.7783120592)
+ 0.10(1)
+ 0.05(0.4935704827)
+ 0.05(1)
= 0.8311053330
~~~

This exactly matches the exported score to displayed precision. It also satisfies every high-confidence rule.

## 11. Cross-validation and leakage audit

### v1

The historical v1 pipeline placed imputation, scaling, and modelling inside scikit-learn pipelines fit on each training fold. Its strict univariate selector was also fit inside each training fold. Those are appropriate leakage controls.

Tier2 is a pre-defined candidate set. The Tier4 top-k feature sets were constructed before CV using full-data gene statistics. However, the current final scoring script removes all non-Tier2_fixed v1 events, so those Tier4 feature sets do not contribute to the current score despite claims in saved documentation.

### strict v2

Strict-v2 feature selection is performed within each outer training split. Hyperparameter selection is nested inside the training data and uses balanced accuracy. This design appropriately separates the outer test fold from selection and tuning.

### Remaining interpretation cautions

- The score is a prioritisation composite, not an unbiased probability of biological truth.
- Repeated folds, nested top-k lists, and models are correlated pieces of evidence, not independent replications.
- Gene-level lmm evidence is derived from the same study and therefore is not external validation.
- No cross-validated ROC, balanced accuracy, or other predictive-performance metric enters the score. Performance tables are read by the current script but assigned to unused variables.

## 12. Double-counting and correlated-evidence audit

Several weighted components reuse related underlying events:

1. Combined selection frequency rewards how often a gene is selected.
2. Model-support count rewards breadth across model families using those same selected events.
3. Sparse support again rewards selection by L1 or elastic net.
4. Importance-rank score summarises the ranks of the same selected rows.
5. Sign consistency summarises the values of the same selected rows.
6. v2 support is triggered by any positive v2 selection frequency, which is already part of combined frequency.

This is intentional multi-criterion prioritisation, but the components are not statistically independent. The score should not be presented as accumulating seven independent confirmations.

The strongest structural duplication is in strict v2: nested top-k sets repeatedly count the same high-ranked gene, and the majority baseline can add selections without using expression features for prediction.

## 13. Comparison with existing method descriptions

| Existing claim or implication | Code-audit finding | Assessment |
|---|---|---|
| Final evidence integrates Tier2_fixed, Tier4 top-k, and strict v2 | Current final script filters v1 to Tier2_fixed; v2 remains | Material discrepancy |
| Tier4 top-k supports the final score | Tier4 rows are removed before final v1 aggregation | Incorrect for current exported score |
| Sign consistency is agreement across signed linear coefficients | Random-forest positive importances are also counted | Material semantic discrepancy |
| v2 supplies stability across nested selection | True, but denominator includes all method/k/model/fold events and the majority baseline | Needs fuller definition |
| Five model families provide model support | Correct; baseline is excluded from model-family count | Correct with baseline caveat |
| Score is based on predictive performance | ROC/performance files do not enter the score | Would be incorrect |
| High/moderate cutoffs are 0.70/0.45 with frequency and support gates | Matches code | Correct |

The run manifest accurately records many source paths and parameters, but its source-strategy description is misleading for the final script because it lists Tier4 and strict-v1 evidence that are filtered out.

## 14. Implementation bugs versus methodological choices

### Likely implementation or documentation defects

- Random-forest importance is allowed into a signed-direction consistency metric.
- Majority-baseline events can increase v2 gene-selection frequency.
- Saved documentation says Tier4/v1-strict evidence contributes, while the final script filters it out.
- Raw event-level inputs needed for full provenance are absent from the current tree.
- The final export omits several derived score components, making independent reconstruction harder than necessary.

### Defensible methodological choices that require explicit wording

- Using a region-global denominator rather than a gene-eligibility denominator.
- Treating nested top-k recurrence as stability evidence.
- Selecting the top 20% for dense linear and forest models.
- Giving v1 and v2 source-level mean ranks equal weight.
- Region-wise scaling of model support, ranks, and lmm evidence.
- Composite weights of 40%, 15%, 10%, 15%, 10%, 5%, and 5%.
- Confidence cutoffs and their additional frequency/support gates.

These are design choices rather than conventional inferential guarantees. They should be labelled as pre-specified prioritisation rules.

## 15. Recommended corrections before reporting or future reruns

### Documentation-only corrections for the current frozen results

1. State that the exported v1 contribution is **Tier2_fixed only**.
2. Define the executed denominators as 125 v1 and 180 v2 region-level events.
3. Explain that v2 frequency spans the full nested method/top-k/model/fold grid and is not conditional on gene eligibility.
4. Describe the score as a heuristic, correlated-evidence prioritisation index, not a probability or significance statistic.
5. State that classification performance does not enter the gene score.
6. Avoid calling the current sign metric “signed linear-model agreement.”
7. Retain the statistical field name lmm_evidence and document its upstream derivation.

### Code changes for a future version, not retroactive changes to this result

1. Exclude random forest and other unsigned models from direction/sign consistency.
2. Exclude the majority baseline from feature-level selection evidence.
3. Decide explicitly whether Tier4 and strict-v1 sources should contribute; make code, manifest, and method text agree.
4. Export event counts and every normalised component used in the weighted sum.
5. Save immutable fold-level input tables with checksums beside each final ranking.
6. Consider reporting two v2 stability measures: unconditional grid support and conditional-on-eligibility support.
7. Add unit tests that reconstruct at least one real gene and verify denominators, component scaling, final score, and confidence label.

## Final audit conclusion

The saved score values are internally reproducible from the current implementation, and the worked CG25 example matches exactly. The audit does **not** support describing the score as evidence integrated across all v1 Tier4 and strict feature sets, because the current final code retains only v1 Tier2_fixed plus strict v2. It also does not support interpreting the current sign-consistency field as purely signed cross-model agreement. The ranking remains usable as an exploratory, within-study gene-prioritisation score provided these implementation-specific definitions and limitations are made explicit.
