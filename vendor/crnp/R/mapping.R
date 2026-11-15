#' Audit feature-to-protein mappings
#'
#' Applies a conservative one-feature-to-one-frozen-protein contract to a
#' resolved mapping table. A feature assigned to multiple frozen proteins is
#' excluded. Multiple eligible features assigned to one protein are retained
#' for downstream aggregation.
#'
#' @param mapping Data frame containing feature and protein identifiers.
#' @param frozen_proteins Character vector of frozen module proteins.
#' @param feature_col Name of the feature identifier column.
#' @param protein_col Name of the resolved protein identifier column.
#' @param tier_col Optional name of a mapping-tier column.
#' @param allowed_tiers Optional accepted mapping tiers.
#'
#' @return Object of class `crnp_mapping_audit`.
#' @export
audit_mapping <- function(mapping, frozen_proteins,
                          feature_col = "feature_id", protein_col = "protein_id",
                          tier_col = NULL, allowed_tiers = NULL) {
  mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)
  required <- c(feature_col, protein_col)
  if (!all(required %in% names(mapping))) {
    stop("mapping lacks required feature or protein columns", call. = FALSE)
  }
  m <- unique(mapping)
  m$feature_id <- as.character(m[[feature_col]])
  m$protein_id <- as.character(m[[protein_col]])
  m$in_frozen_module <- !is.na(m$protein_id) & m$protein_id %in% frozen_proteins
  if (!is.null(tier_col)) {
    if (!tier_col %in% names(m)) stop("tier_col is absent from mapping", call. = FALSE)
    m$mapping_tier <- as.character(m[[tier_col]])
    m$tier_allowed <- if (is.null(allowed_tiers)) TRUE else
      m$mapping_tier %in% allowed_tiers
  } else {
    m$mapping_tier <- NA_character_
    m$tier_allowed <- TRUE
  }
  counts <- tapply(m$protein_id[m$in_frozen_module],
                   m$feature_id[m$in_frozen_module],
                   function(x) length(unique(x)))
  ambiguous <- names(counts)[counts > 1L]
  m$ambiguous_feature <- m$feature_id %in% ambiguous
  m$eligible <- m$in_frozen_module & m$tier_allowed & !m$ambiguous_feature
  eligible <- unique(m[m$eligible, c("feature_id", "protein_id", "mapping_tier"),
                       drop = FALSE])
  protein_counts <- table(factor(eligible$protein_id, levels = frozen_proteins))
  protein_summary <- data.frame(
    protein_id = frozen_proteins,
    eligible_feature_count = as.integer(protein_counts),
    structurally_measurable = as.integer(protein_counts) > 0L,
    stringsAsFactors = FALSE
  )
  summary <- data.frame(
    frozen_proteins = length(frozen_proteins),
    mapped_proteins = sum(protein_summary$structurally_measurable),
    structural_coverage = mean(protein_summary$structurally_measurable),
    input_features = length(unique(m$feature_id)),
    eligible_features = length(unique(eligible$feature_id)),
    ambiguous_features = length(unique(m$feature_id[m$ambiguous_feature])),
    stringsAsFactors = FALSE
  )
  structure(list(
    feature_manifest = m,
    eligible_mapping = eligible,
    protein_summary = protein_summary,
    summary = summary
  ), class = "crnp_mapping_audit")
}

#' Print a protein-mapping audit
#' @param x A mapping audit.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.crnp_mapping_audit <- function(x, ...) {
  s <- x$summary[1L, ]
  cat("CRNP mapping audit\n")
  cat("  Frozen proteins:", s$frozen_proteins, "\n")
  cat("  Structurally measurable:", s$mapped_proteins, "\n")
  cat("  Structural coverage:", sprintf("%.1f%%", 100 * s$structural_coverage), "\n")
  cat("  Excluded ambiguous features:", s$ambiguous_features, "\n")
  invisible(x)
}

#' Aggregate assay features to one value per protein
#'
#' Distinguishes participant-level missing measurements from proteins absent
#' from an assay. Eligible features are filtered by participant completeness,
#' optionally median-imputed, standardized separately, and combined within
#' protein by the participant-wise median or mean.
#'
#' @param feature_matrix Numeric participants by assay-feature matrix.
#' @param mapping Resolved feature-to-protein mapping data frame.
#' @param frozen_proteins Frozen module proteins.
#' @param feature_col Mapping feature column.
#' @param protein_col Mapping protein column.
#' @param tier_col Optional mapping-tier column.
#' @param allowed_tiers Optional accepted mapping tiers.
#' @param min_finite_fraction Minimum finite participant fraction for a feature.
#' @param sample_missing `"median"` or `"error"`.
#' @param aggregate `"median"` or `"mean"` for multiple features per protein.
#'
#' @return List containing the protein matrix and complete audit tables.
#' @export
aggregate_features <- function(feature_matrix, mapping, frozen_proteins,
                               feature_col = "feature_id",
                               protein_col = "protein_id", tier_col = NULL,
                               allowed_tiers = NULL,
                               min_finite_fraction = 0.8,
                               sample_missing = c("median", "error"),
                               aggregate = c("median", "mean")) {
  sample_missing <- match.arg(sample_missing)
  aggregate <- match.arg(aggregate)
  x <- .crnp_validate_matrix(feature_matrix, "feature_matrix", require_finite = FALSE)
  if (!is.finite(min_finite_fraction) || min_finite_fraction <= 0 ||
      min_finite_fraction > 1) {
    stop("min_finite_fraction must be in (0, 1]", call. = FALSE)
  }
  audit <- audit_mapping(
    mapping, frozen_proteins, feature_col, protein_col, tier_col, allowed_tiers
  )
  eligible <- audit$eligible_mapping
  eligible <- eligible[eligible$feature_id %in% colnames(x), , drop = FALSE]
  if (!nrow(eligible)) stop("No eligible mapped assay features", call. = FALSE)
  feature_names <- unique(eligible$feature_id)
  finite_fraction <- vapply(feature_names, function(f) {
    mean(is.finite(x[, f]))
  }, numeric(1))
  retained <- feature_names[finite_fraction >= min_finite_fraction]
  if (!length(retained)) stop("No features meet the completeness threshold", call. = FALSE)
  x2 <- x[, retained, drop = FALSE]
  imputed_count <- stats::setNames(integer(ncol(x2)), colnames(x2))
  for (j in seq_len(ncol(x2))) {
    bad <- !is.finite(x2[, j])
    imputed_count[[j]] <- sum(bad)
    if (any(bad)) {
      if (sample_missing == "error") {
        stop("Participant-level missing values remain in an eligible feature",
             call. = FALSE)
      }
      x2[bad, j] <- stats::median(x2[!bad, j], na.rm = TRUE)
    }
  }
  feature_sd <- apply(x2, 2L, stats::sd)
  nonconstant <- is.finite(feature_sd) & feature_sd > 0
  x2 <- x2[, nonconstant, drop = FALSE]
  if (!ncol(x2)) stop("All eligible features have zero variance", call. = FALSE)
  feature_z <- .crnp_standardize(x2)$matrix
  eligible <- eligible[eligible$feature_id %in% colnames(feature_z), , drop = FALSE]
  proteins <- frozen_proteins[frozen_proteins %in% eligible$protein_id]
  protein_matrix <- vapply(proteins, function(protein) {
    features <- unique(eligible$feature_id[eligible$protein_id == protein])
    values <- feature_z[, features, drop = FALSE]
    if (ncol(values) == 1L) return(values[, 1L])
    if (aggregate == "median") {
      apply(values, 1L, stats::median)
    } else {
      rowMeans(values)
    }
  }, numeric(nrow(feature_z)))
  if (is.null(dim(protein_matrix))) {
    protein_matrix <- matrix(protein_matrix, ncol = 1L,
                             dimnames = list(rownames(feature_z), proteins))
  }
  protein_matrix <- .crnp_standardize(protein_matrix)$matrix
  rownames(protein_matrix) <- rownames(x)
  feature_qc <- data.frame(
    feature_id = feature_names,
    finite_fraction = unname(finite_fraction[feature_names]),
    retained_by_completeness = feature_names %in% retained,
    retained_nonconstant = feature_names %in% colnames(feature_z),
    imputed_values = ifelse(
      feature_names %in% names(imputed_count),
      imputed_count[match(feature_names, names(imputed_count))], 0L
    ),
    stringsAsFactors = FALSE
  )
  list(
    protein_matrix = protein_matrix,
    observed_proteins = colnames(protein_matrix),
    structural_missing_proteins = setdiff(frozen_proteins, colnames(protein_matrix)),
    mapping_audit = audit,
    feature_qc = feature_qc,
    aggregation = aggregate,
    min_finite_fraction = min_finite_fraction,
    sample_missing = sample_missing
  )
}
