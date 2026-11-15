# Environment installation

The tested numerical environment is Python **3.11.7**, NumPy **1.26.4**, pandas **2.1.4**, SciPy **1.17.0**, scikit-learn **1.2.2**, R **4.4.3**, and WGCNA **1.73**. These values come from the frozen manuscript environment. Numerical BLAS/OpenMP pools are restricted to one thread in the runners.

Use an isolated Python environment and a repository checkout:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements-manuscript-lock.txt
python -m pip install --no-deps -e .
```

For demonstration-only use, `requirements.txt` permits a broader compatible environment; its results are not asserted to be bitwise manuscript reruns. Use Python 3.11 with the supplied scikit-learn 1.2.2 constraint. The lock file fixes the five numerical environment versions; plotting/helper versions are explicitly separated in that file. This is not a complete OS/BLAS/Bioconductor snapshot or a claim that hardware-dependent timings will match.

**Editable installation is required** when installing the CLI: resources remain in the checkout next to `src/`. The scripts also work directly from the checkout without installing the CLI once dependencies are present. A standalone wheel without the repository's `metadata/`, `R/` and `vendor/` trees is not a supported distribution.

For R, install R 4.4.3 for strict reruns and configure its executable on PATH, or use the `RSCRIPT` environment variable. The runner sources CRNP v0.2.0 directly; a separately installed `crnp` package is unnecessary.

```powershell
Rscript R/install_dependencies.R
```

The helper installs WGCNA 1.73, data.table 1.18.2.1, jsonlite 2.0.0 and microbenchmark 1.5.0, using remotes for archived versions if needed, plus WGCNA's Bioconductor dependencies. It needs internet access, a writable R library, and possibly the R build toolchain for the selected operating system. Use a new user library if you do not want to alter an existing research environment. Bioconductor's complete dependency tree is not frozen here. If installation fails, its ordinary R error identifies the missing dependency; the analysis never silently substitutes another CRNP implementation.

An existing R environment may be selected with `R_LIBS`/`R_LIBS_USER`. These variables belong to your local shell; do not commit machine-specific configuration. `LC_ALL=C` is set only in the worker subprocess; on Windows the worker explicitly selects a UTF-8 character locale. Filenames are passed through UTF-8 JSON. Chinese-path operation is covered in local validation.

Run the following from the repository root:

```powershell
python -B -m unittest discover -s tests -v
python -B scripts/check_r_contracts.py
python -B examples/run_demo.py --bootstrap 50
python -B scripts/redraw_paper_figures.py
```

Pure Python tests and aggregate plotting do not require R. The full demonstration and R contracts do. Demonstration runs allow and record environment drift; governed-data runs enforce the five Python numerical versions and R/WGCNA versions by default. `--allow-version-drift` is an explicit override and is reported in run status. That switch does not validate new versions automatically.

The local release validation used existing isolated environments. The dependency installers have not been tested by rebuilding every dependency on a clean operating system. See `validation/VALIDATION_SUMMARY.json` for the exact observed checks and limitations.
