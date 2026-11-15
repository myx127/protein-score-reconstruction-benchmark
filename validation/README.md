# Local validation evidence

All quantitative validation data generated during release preparation were **independent synthetic values and identifiers**. No HMACS/Emory participant matrix or real-data fitted model is included in these reports. The release task did not refit the actual cohorts to establish complete manuscript reproduction.

| Evidence file | Meaning |
| --- | --- |
| `VALIDATION_SUMMARY.json` | Consolidated actual checks, counts, environments and limits. |
| `python_tests.txt` | Passing scientific/data-boundary tests, including aggregation/paired-bootstrap contracts. |
| `frozen_driver_parity.json` | Comparison with the exact original aggregate/bootstrap driver on a 72-setting synthetic prediction fixture. |
| `paired_assay_point_parity.json` | Independent Python rank-correlation/Fisher-z recomputation of the R synthetic endpoint. |
| `R_CONTRACTS.json` | Five CRNP projection contracts plus parity against an installed frozen CRNP 0.2.0 when available. |
| `synthetic_recovery_status.json` | Two directions, 24 settings, four methods and 50 conditional bootstrap draws. |
| `synthetic_emory_status.json` | Six paired branches, 35 generated participants, 18 cells and 50 draws. |
| `synthetic_timing_status.json` | 24 workloads, 72 warm-ups, 360 formal records, five repetitions, interleaved languages. |
| `editable_install_check.json` | Isolated editable build, resource lookup and CLI help; existing research environment unmodified. |
| `input_encoding_check.json` | R input preserves leading zeros in UTF-8-BOM compressed CSV. |
| `public_release_audit.json` | Public CSV schemas, file counts, text scan and allowed-content inventory. |

Numerical parity does not establish independent biological validation, and synthetic execution does not establish the manuscript's effect sizes. The broader interpretation limits are in `docs/REPRODUCIBILITY.md`. Local run logs/participant-like toy matrices remain outside this public package; only summaries are copied here.
