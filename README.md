# Protein network score reconstruction benchmark

Code accompanying **Benchmarking Protein Network Score Reconstruction under Structured Feature Missingness**. Prepared locally for public release, version **1.0.0**; this is not a claim of publication or conference acceptance. The CRNP algorithm included here remains version **0.2.0**.

[中文使用说明](README_中文.md) · [Installation](docs/INSTALL.md) · [Input contracts](docs/INPUT_SCHEMAS.md) · [Data access](docs/DATA_ACCESS.md) · [Reproducibility](docs/REPRODUCIBILITY.md)

The public package contains frozen numerical functions, portable Python/R runners, independent synthetic examples, protein/mask metadata, and manuscript aggregate results. **No HMACS or Emory participant-level quantitative matrices, participant scores, or fitted patient-data models are distributed.**

## Quick start

Use a checkout of this repository. Python 3.11 and R are needed; install the R dependencies as described in [INSTALL.md](docs/INSTALL.md). From the repository root:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python -m pip install --no-deps -e .
Rscript R/install_dependencies.R
python -B -m unittest discover -s tests -v
python -B examples/run_demo.py --bootstrap 50
```

If `Rscript` is not on PATH, set `$env:RSCRIPT` to its executable path. Use `R_LIBS` to select an existing R package library when necessary. The demonstration writes only generated synthetic data to `outputs/synthetic_demo/`, including 24 recovery settings and a six-branch, 35-person synthetic paired-assay example. Its three toy modules and effect sizes are independent of the manuscript. Use a fresh `--output` directory for every run.

```powershell
python -B examples/run_demo.py --output outputs/demo_with_timing --bootstrap 50 --with-timing
python -B scripts/check_r_contracts.py
python -B scripts/redraw_paper_figures.py
python -B tools/verify_release.py
```

Timing is optional and runs 24 workloads × 3 methods × 5 repetitions, with 72 additional warm-ups. It preserves a persistent R process and randomized interleaving across languages. Machine-specific runtimes will differ from the manuscript.

## What can be reproduced

| Access level | Supported work |
| --- | --- |
| Public files only | Run numerical tests, synthetic recovery/paired-assay/timing demos, check aggregate results, redraw figures from published aggregate tables. |
| Approved frozen HMACS matrices | Refit the 72-setting, two-direction recovery benchmark and 2,000-draw conditional participant bootstrap; rerun the 24-workload timing experiment. |
| Approved HMACS definition matrix plus permitted six-branch Emory inputs | Refit the paired-assay scores and 10,000-draw concordance bootstrap, conditional on the documented mapping/preprocessing. |

The analysis starts from the **frozen, preprocessed representation**. It does not reconstruct raw mass-spectrometry processing, original cohort recruitment, or joint module discovery. The two standardized HMACS batches do not substitute for the separate source-definition matrix used in the Emory branch. See [REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for the exact scope and local validation evidence.

## Governed-data analysis

Keep quantitative inputs, local configurations, predictions, and fitted objects **outside this repository**. Obtain the approved files before running; the example configurations are templates, not downloaders. Adjust local paths in an external copy of each example:

```powershell
python -B scripts/run_benchmark.py recovery --config ../governed_inputs/hmacs.local.json --private-workdir ../private_runs/recovery_run1
python -B scripts/run_benchmark.py timing --config ../governed_inputs/hmacs.local.json --private-workdir ../private_runs/timing_run1 --repetitions 5
python -B scripts/run_benchmark.py emory --config ../governed_inputs/emory.local.json --private-workdir ../private_runs/emory_run1
```

Examples are in `config/`. Defaults use the public frozen weights, masks, directions and seeds. HMACS matrix checksums are enforced by default. `--allow-version-drift` records an explicit environment deviation; it does not establish manuscript reproduction. Disabling input fingerprint checks is a new-data analysis, not an exact rerun.

## Repository contents

```text
src/psbench/          Frozen numerical functions and portable runners
vendor/crnp/         Preserved CRNP v0.2.0 R source, help and license
R/                   R worker and dependency installer
examples/            Independent synthetic-data demonstration
scripts/             Analysis, R checks and aggregate figure entry points
config/              Governed input examples; no real paths or data
metadata/            Protein weights, feature mappings, masks, seeds, hashes
results/published/   Frozen manuscript aggregate tables
paper_figures/       Manuscript figure PDF snapshots
tests/               Scientific and data-boundary checks
validation/          Sanitized evidence from local release validation
docs/                Methods, data access, reproducibility, upload guide
```

The benchmark uses normalized fixed-weight scores. Source residual ridge and CRNP are source-trained methods. Within-CF ridge uses target training-fold full scores and is a **diagnostic with stronger information**, not a deployable comparator or a mathematical bound. HMACS uses raw direct scoring; Emory uses source-calibrated, branch-standardized direct scoring. The Emory intervals include zero and do not establish superiority or equivalence.

## Citation, license and contact

See [CITATION.cff](CITATION.cff). Authors: Xinyue Wang, Yuxi Mao and Yan Zeng. A repository URL/archival DOI can be added after the owner creates the public repository and release; no identifier is invented here. Code is MIT licensed; see [LICENSE](LICENSE) and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). The code license grants no participant-data access rights.

HMACS access requests: **Yan Zeng, zengyan68@wust.edu.cn**, subject to the applicable cohort/institutional review and data-use conditions. Requests do not guarantee approval. Third-party data remain subject to their original conditions.
