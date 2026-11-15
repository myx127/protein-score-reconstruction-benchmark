# v1.0.0 - Reproducible benchmark code

Initial public software release accompanying *Benchmarking Protein Network Score Reconstruction under Structured Feature Missingness*. The benchmark package version is 1.0.0; the preserved CRNP algorithm remains version 0.2.0. Software publication does not imply manuscript publication or conference acceptance.

- Frozen Python numerical functions and original CRNP R source.
- Independent synthetic recovery, paired-assay and interleaved timing demonstrations.
- Configurable governed-input runners, public feature/coverage metadata and manuscript aggregate results.
- Environment instructions, input contracts, data access statement, MIT license, citation metadata and SHA256 release manifest.

Validation: 15 Python checks; five R projection contracts plus installed-CRNP parity; a 24-setting synthetic recovery example and 35-person six-branch synthetic paired-assay example with 50 bootstrap draws; 24 timing workloads, 72 warm-ups and 360 formal records; aggregate figure redrawing and Unicode-path execution. The portable aggregate/bootstrap adapter agrees with the original frozen driver to floating-point tolerance on a 72-setting synthetic fixture.

No HMACS or Emory participant-level matrices, participant scores or real-data fitted models are included. HMACS access requests go to the corresponding author, Yan Zeng (zengyan68@wust.edu.cn), subject to applicable cohort/institutional approval and data-use conditions. Third-party data retain their original terms. This release preparation did not refit actual cohort matrices to claim complete manuscript reproduction.

Start with `README.md` or `README_中文.md`. Verify the package with `python -B tools/verify_release.py`. For generated synthetic data, follow `docs/INSTALL.md` and run `python -B examples/run_demo.py --bootstrap 50`.
