# Numerical and design contracts

Let a module have frozen signed weights `w` and protein values `X`. The full reference is `F = X_M w_M / sum(abs(w_M))`. For the observed set `O`, the raw direct score is `D = X_O w_O / sum(abs(w_O))`. Their difference includes both missing weighted contributions and a changed normalization denominator. The reference is a fixed network score, not a measured clinical truth or a newly estimated eigengene.

## Source residual ridge

Regress source `F-D` on the observed source protein values, with an unpenalized intercept. The eight alpha values are `0.0001,0.001,0.01,0.1,1,10,100,1000`. Cross-validation is shuffled, with `min(5,max(2,n//20))` folds and the registered deterministic seed. Select the **largest** alpha whose mean validation MSE is within one standard error of the minimum; use the original `1e-15` numerical tolerance. Refit on all complete source rows with scikit-learn's Cholesky ridge solver. Prediction is target `D + predicted residual`, without target full-score calibration. The R package's internal `ridge` alternative is not used in place of this method.

## Within-target cross-fitting diagnostic

Use five outer target folds; tune ridge within each outer training fold using the same selector and a derived seed. Target full reference scores are available in those training folds. Each evaluated participant is held out from its fitted model. This is a stronger-information diagnostic for remaining recoverable discrepancy. Its empirical error gap can be negative in a setting; it is not a theoretical upper bound and need not be achieved by source-only models. No deployment or prospective prediction claim follows.

## CRNP v0.2.0

The included R source is preserved byte for byte. It fits source biweight midcorrelations, regularized missing-score reconstruction, source calibration and the coverage attenuation `a(c)=max(0,min(1,c*(c-0.20)))`. Its six regularization values are `0.0001,0.001,0.01,0.1,1,10`, with three-fold source calibration and eight masks: random, hub-depleted, hub-enriched and topology-blocked, twice each. The original full fitting API also constructs auxiliary method models; these remain included in fitting and in CRNP timing.

For HMACS recovery, input batches are already frozen standardized matrices. Projection explicitly sets observed-feature centers to zero and scales to one, with `standardize_score=FALSE`. Because the CRNP score is calibrated on the standardized source full score, it is transformed back with the **source full-score mean and sample SD**. Target full scores are used only for evaluation, not calibration. Correcting this scale is essential for meaningful normalized MSE.

For the Emory evaluation, the preserved `aggregate_features` and `compare_methods` machinery implements conservative feature mapping, median aggregation and source calibration. Both direct and CRNP scores are standardized by branch. Thus the direct comparator here differs from HMACS's raw direct score. Cross-platform concordance tests agreement of those branch scores, not accuracy against a known complete target score.

## Estimands and uncertainty

For each recovery setting/method, `risk = mean((F-P)^2) / var(F,ddof=1)`. First average 36 settings equally within each source-target direction, then average the two directions equally. Compute gain `R_direct-R_method`, relative gain `(R_direct-R_method)/R_direct`, and the fraction of empirical diagnostic headroom **after** aggregating risks. Do not average setting-level ratios or weight the directions by participant number. Fractions with a nonpositive denominator remain undefined.

Primary intervals use 2,000 registered target participant draws, sharing the draw across all settings and methods within each target batch. All records for a sampled participant travel together. Fitted models, preprocessing and predictions are fixed. The efficient adapter resamples cached squared errors and recomputes the target sample variance, which is algebraically equivalent to resampling the fixed score vectors. It does not capture uncertainty from module discovery, preprocessing or source model fitting. Settings share proteins and participants and are not independent biological replicates.

Emory's endpoint averages, equally across evaluable cells, `atanh(Spearman_CRNP)-atanh(Spearman_direct)`, clipping correlations at +/-0.999999. A cell is module × platform pair × biofluid. The same R `sample.int` draw is used for all six branches among 35 paired participants; scores stay fixed and correlations are recalculated. The frozen bootstrap seed is 2026081601 with 10,000 draws. Overall, CSF and plasma summaries are reported. Intervals crossing zero do not establish equivalence. No formal source-ridge-versus-CRNP superiority test is asserted.

## Timing

The paper uses three modules (40,159,285 proteins), four masks, two directions: 24 workloads. Five formal repetitions per method yield 360 records, with 72 warm-ups. Common data loading, extraction, startup and preprocessing are outside timers. Source score/residual construction, all tuning, calibration, the complete CRNP fitting API and fresh fitting are inside fit timers; score timers evaluate a whole target batch with a fitted model. Direct scoring has a structural zero fit cost.

Python and the resident R worker receive a common randomized interleaved command order; numeric pools are single-threaded. Warm-ups freeze power-of-two prediction block sizes targeting about 0.1 seconds, bounded from 32 to 131072 calls. Each timed result is checked against an untimed fit on identical inputs. Repetitions are summarized by workload median; the overall median/IQR is across those workload medians. The IQR represents workload heterogeneity, not a confidence interval. Implementation, API and language overhead remain part of this practical cost comparison; it does not isolate an algorithmic complexity theorem.

## Representation and interpretation limits

The HMACS batches are participant-disjoint, but the upstream module representation/preprocessing was frozen from a joint analysis. Recovery validates transfer conditional on that representation; it is not an independently discovered full pipeline or clinical outcome prediction study. MGH and Sissala provide structural protein-set masks, not quantitative values in the HMACS recovery experiment. Paired-assay Emory analysis provides a separate, limited concordance evaluation with no demonstrated improvement. Synthetic outputs validate software execution and contracts; manuscript aggregate tables remain the frozen scientific evidence unless governed-data reruns are separately completed and reconciled.
