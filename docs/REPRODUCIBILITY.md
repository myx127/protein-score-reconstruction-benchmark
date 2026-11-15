# Reproducibility scope and local evidence

The release supports public execution **without participant-data access**, and governed-data reruns **after approved inputs are provided**. These are different levels of evidence.

| Level | Meaning | Included/tested |
| --- | --- | --- |
| Public software execution | All three methods and the stronger-information diagnostic run on independent generated values. | Synthetic two-batch recovery, paired six-branch concordance, bootstrap, and full 24-workload timing. |
| Frozen aggregate evidence | Paper summary numbers and settings can be checked and plotted without data access. | Six aggregate CSV files, original figure PDFs, portable redrawing script, reconciliation tests. |
| Exact governed-input rerun | Refit/evaluate from approved frozen participant matrices in the frozen environment. | Configurable runners and file-hash guards provided; **not performed on real data during this release-preparation task**. |
| Full upstream clinical/raw-data reproduction | Recruitment, raw DIA processing, normalization and joint network discovery independently reproduced. | Outside this package's scope. Frozen representations are inputs. |

## Preserved implementation

`src/psbench/frozen.py` retains 17 original function/class bodies, including the ridge selector, nested cross-fitting, metrics, prediction generation, participant resampling and percentile inference. `metadata/code_provenance.json` records their original body hashes and the byte hashes of all six CRNP R source files. Tests verify these hashes. Public metadata preserve source weight strings to avoid additional numeric rounding during export; the numeric round-trip was checked against the original metadata.

The original aggregate/bootstrap driver hard-coded 36 settings per direction. The portable adapter accepts the smaller synthetic design and computes the same equal-cell/equal-direction estimand. Its cached-squared-error bootstrap was compared against the **exact original driver functions**, with independently simulated predictions in the complete 72-setting layout and five paired draws. See `validation/frozen_driver_parity.json` for observed differences. This validates the adapter's numerical equivalence for finite locked predictions; it does not assert a real-cohort refit.

CRNP projection checks cover the effective-weight/calibration identity, invariance to unobserved target feature changes, protein-column reordering, attenuation, and inversion with source full-score moments. When CRNP 0.2.0 is already installed, the R check also compares predictions against that installed implementation. No installed package is required for normal execution because the original source is vendored.

## Validation evidence

`validation/VALIDATION_SUMMARY.json` records the checks actually completed locally, the observed environment, synthetic counts and the release boundary. Local validation covered the Python scientific tests, R numerical contracts, complete recovery/paired-assay demonstration with 50 bootstrap draws, Unicode file paths, the full 72-warm-up/360-record timing demonstration, aggregate figure redrawing, and release checks. Each generated participant-like value in those demonstrations was synthetic. Real HMACS/Emory quantitative matrices were neither copied into the release nor refitted to claim full reproduction.

The small demonstration draw count is intentionally for software validation. It is not suitable for a precise inferential claim. Governed configurations default to 2,000 recovery and 10,000 paired-assay draws, as in the paper. Timing summaries remain environment-dependent. The source numerical environment and plotting environment are recorded separately where applicable; no clean-OS dependency rebuild or cross-platform validation is claimed.

The package's dependency installer supplies version targets, not a complete immutable OS/container image. Strict numerical checks cover Python, NumPy, pandas, SciPy, scikit-learn and R/WGCNA. R dependency details are written in each run's `R_environment.json`. A recorded environment override is not evidence that another environment reproduces the manuscript.

## Manuscript aggregate sanity checks

HMACS includes 72 related settings (9 modules × 4 masks × 2 directions), with 60 improvements and 12 deteriorations for each source-trained method. The public full-score reference error is about 0.047856 for direct, 0.043384 for source residual ridge, 0.044694 for CRNP, and 0.035118 for Within-CF. These are descriptive checks of the frozen tables, not new measurements. The paired-assay mean difference is approximately -0.002266 and its interval includes zero. No formal ridge-over-CRNP superiority or paired-assay equivalence conclusion is encoded.

`paper_figures/` contains the manuscript PDF snapshots. Redrawn figures use the same numerical aggregate inputs but a portable layout/open font and are labeled `_redrawn`; identical typography or panel geometry is not asserted. The workflow drawing is regenerated from the documented design rather than inferred from participant data.

## After a governed rerun

Keep all outputs external. First inspect failure files, settings/sample counts, frozen fingerprints, environment reports, fold-overlap audits and finite prediction coverage. Compare `aggregate/point_estimates.csv` with the included tables using numeric tolerances and the same input versions; distinguish score-level floating-point tolerance from inferential differences. Check bootstrap counts and scope separately. Record discrepancies before describing an exact rerun. Selecting new weights, masks, normalization, participant samples, seeds or target calibration defines a different analysis and should receive a separate configuration/version.
