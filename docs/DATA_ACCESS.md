# Participant data access and public release boundary

HMACS denotes the Hubei Memory and Aging Cohort Study. **HMACS participant-level data are not publicly distributed in this code repository.** Interested researchers may contact the corresponding author, **Yan Zeng (曾燕), zengyan68@wust.edu.cn**, describing the research purpose, requested variables/files, institution, authorized users, and intended storage/analysis arrangements. Access is subject to the applicable cohort/institutional review and data-use conditions; a request does not guarantee approval.

For an exact recovery rerun, request the approved frozen matrices named in `metadata/expected_governed_inputs.json`, together with their preprocessing description. The recovery benchmark has two participant-disjoint batches of 722 and 833 participants and 1,222 protein columns. Matrix hashes identify the exact frozen compressed files; different compression or line endings may alter a file hash even when numeric values agree. Do not disable a failed hash check and describe the run as exact without independently establishing equivalence.

For the Emory branch, also request the separate governed HMACS **source-definition matrix** used with the raw UniProt module weights. Do not reconstruct it by concatenating the two standardized recovery batches. Emory participant inputs must be obtained and used under their original provider's terms. This repository does not host them or grant their redistribution. Its assay/protein mapping records are feature metadata only.

The public release includes group-level module weights, protein membership and assay mapping, deterministic seed tables, and frozen **aggregate** result tables. No real participant rows, cross-assay sample links, per-person scores, fitted source models/covariance objects, clinical phenotypes, request credentials or data-access tokens are included. Synthetic examples use independent values and generated identifiers.

Place governed matrices, local configuration, scores, models, logs, bootstrap intermediates and any joined clinical information **outside the repository root**. The governed-data entry points reject input/output paths within the checkout. This boundary prevents routine accidental mixing; it does not determine whether an arbitrary file is safe to share. The `.gitignore` file and release checker provide additional mechanical checks and do not override provider/institutional requirements.

Suggested manuscript statements, after filling in the actual published code URL:

> Analysis code, independent synthetic examples, frozen protein/coverage metadata, and aggregate results are available at [insert the actual repository and release identifier]. HMACS participant-level data are not publicly available and may be requested from the corresponding author, Yan Zeng (zengyan68@wust.edu.cn), subject to applicable cohort and institutional approvals and data-use conditions. Third-party participant data remain subject to their original access conditions.

Until a repository/release is actually public, use a local preparation statement in drafts; do not invent a URL, DOI or unrestricted data availability.
