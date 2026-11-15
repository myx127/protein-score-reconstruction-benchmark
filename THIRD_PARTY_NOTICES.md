# Source and dependency notices

`vendor/crnp/` preserves CRNP version 0.2.0, authored by Yuxi Mao. Its original MIT license and attribution remain in `vendor/crnp/LICENSE` and `vendor/crnp/LICENSE.md`; all R function source files are preserved byte for byte. File hashes are recorded in `metadata/code_provenance.json`.

`src/psbench/frozen.py` contains unchanged function/class bodies extracted from the authors' frozen Stage24 benchmark implementation. A new adapter supplies portable paths, small-design aggregation, efficient conditional bootstrap, and persistent-worker timing. Original driver code that contained private paths is not distributed. Numerical-body hashes are recorded and tested.

NumPy, pandas, SciPy, scikit-learn, threadpoolctl, matplotlib, psutil, R, WGCNA, data.table, jsonlite, microbenchmark, and their dependencies are installed separately and retain their own licenses. They are not relicensed by this repository. Python and R installers are environment helpers, not vendored dependency distributions.

Protein identifiers, module weights, mask membership, assay-to-protein annotations and aggregate result tables are scientific metadata/results for this benchmark. Their release provides no license or permission to obtain or redistribute governed cohort or third-party participant-level data. Original third-party data terms continue to apply.

`paper_figures/` contains the authors' manuscript figure snapshots. `scripts/redraw_paper_figures.py` uses the included aggregate tables to redraw the numerical evidence with a portable open-font layout. No font binaries, source manuscript archives, third-party article PDFs, or participant-level matrices are included.
