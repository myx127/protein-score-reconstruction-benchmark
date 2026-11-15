#' Fit a coverage-robust network projection model
#'
#' Fits CRNP from a complete source protein matrix, frozen signed module
#' weights, and the proteins structurally measurable in a target assay. No
#' target values or phenotypes are used for fitting or tuning.
#'
#' @param source_matrix Numeric source samples by protein matrix with unique
#'   protein column names and finite values.
#' @param weights Named numeric vector of frozen signed protein weights.
#' @param observed_proteins Names of proteins structurally measurable in the
#'   target assay.
#' @param lambdas Candidate source-only ridge parameters.
#' @param calibration_mechanisms Source masking mechanisms used to select the
#'   regularization parameter.
#' @param calibration_replicates Replicates per calibration mechanism.
#' @param folds Number of source-only folds.
#' @param power Signed-adjacency power used only for topology-blocked masks.
#' @param seed Random seed for source masks.
#'
#' @return An object of class `crnp_model`.
#' @export
fit_crnp <- function(source_matrix, weights, observed_proteins,
                     lambdas = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10),
                     calibration_mechanisms = c(
                       "random", "hub_depleted", "hub_enriched",
                       "topology_blocked"
                     ),
                     calibration_replicates = 2L, folds = 3L,
                     power = 8, seed = 93017L) {
  source_matrix <- .crnp_validate_matrix(source_matrix, "source_matrix")
  if (is.null(names(weights)) || anyNA(names(weights)) || anyDuplicated(names(weights))) {
    stop("weights must have unique, nonmissing protein names", call. = FALSE)
  }
  weights <- as.numeric(weights) |> stats::setNames(names(weights))
  if (any(!is.finite(weights))) stop("weights must be finite", call. = FALSE)
  if (!all(names(weights) %in% colnames(source_matrix))) {
    stop("All weighted proteins must be present in source_matrix", call. = FALSE)
  }
  observed_proteins <- names(weights)[names(weights) %in% observed_proteins]
  if (length(observed_proteins) < 2L || length(observed_proteins) >= length(weights)) {
    stop("observed_proteins must contain at least two, but not all, module proteins",
         call. = FALSE)
  }
  if (any(!is.finite(lambdas)) || any(lambdas <= 0)) {
    stop("lambdas must be finite and greater than zero", call. = FALSE)
  }
  ordered_source <- source_matrix[, names(weights), drop = FALSE]
  standardized <- .crnp_standardize(ordered_source)
  source_z <- standardized$matrix
  r_source <- .crnp_robust_cor(source_z)
  observed <- match(observed_proteins, names(weights))
  coverage_fraction <- length(observed) / length(weights)
  recovery_fraction <- .crnp_recovery_fraction(coverage_fraction)

  masks <- list()
  position <- 0L
  for (mechanism_i in seq_along(calibration_mechanisms)) {
    for (replicate in seq_len(calibration_replicates)) {
      position <- position + 1L
      masks[[position]] <- .crnp_make_mask(
        weights, r_source, coverage_fraction,
        calibration_mechanisms[[mechanism_i]],
        seed = seed + mechanism_i * 100L + replicate, power = power
      )
    }
  }
  tuning_crnp <- .crnp_tune_regularization(
    source_z, weights, masks, lambdas, method = "crnp", folds = folds,
    recovery_fraction = recovery_fraction
  )
  tuning_ridge <- .crnp_tune_regularization(
    source_z, weights, masks, lambdas, method = "ridge", folds = folds,
    recovery_fraction = 1
  )
  lambda_crnp <- tuning_crnp$best
  lambda_ridge <- tuning_ridge$best
  truth <- .crnp_z_vector(as.numeric(source_z %*% weights))
  methods <- c("direct_kME", "signed_mean", "subset_PC1", "ridge",
               "CRNP_full", "CRNP", "source_mean_imputed_fixed")
  method_models <- lapply(methods, function(method) {
    beta <- .crnp_method_weights(
      method, r_source, weights, observed, lambda_crnp, lambda_ridge,
      recovery_fraction
    )
    names(beta) <- observed_proteins
    list(
      effective_weights = beta,
      calibration = .crnp_calibration(source_z, truth, beta, observed)
    )
  })
  names(method_models) <- methods
  fidelity <- score_fidelity(
    r_source, weights, observed, lambda_crnp, recovery_fraction
  )

  structure(list(
    version = "0.2.0",
    weights = weights,
    observed_proteins = observed_proteins,
    effective_weights = method_models$CRNP$effective_weights,
    coverage_fraction = coverage_fraction,
    recovery_fraction = recovery_fraction,
    lambda = lambda_crnp,
    lambda_ridge = lambda_ridge,
    source_correlation = r_source,
    source_center = standardized$center,
    source_scale = standardized$scale,
    calibration_intercept = method_models$CRNP$calibration[["intercept"]],
    calibration_slope = method_models$CRNP$calibration[["slope"]],
    correlation_fidelity = fidelity[["correlation_fidelity"]],
    mse_fidelity = fidelity[["mse_fidelity"]],
    tuning_summary = tuning_crnp$summary,
    ridge_tuning_summary = tuning_ridge$summary,
    method_models = method_models
  ), class = "crnp_model")
}

#' Project a fitted CRNP model
#'
#' @param model A fitted `crnp_model`.
#' @param target_matrix Numeric target samples by measured-protein matrix.
#' @param center Optional target reference means.
#' @param scale Optional target reference standard deviations.
#' @param standardize_score Whether to standardize the final score.
#'
#' @return Numeric participant-level CRNP score with target reference
#'   parameters stored as attributes.
#' @export
project_crnp <- function(model, target_matrix, center = NULL, scale = NULL,
                         standardize_score = FALSE) {
  if (!inherits(model, "crnp_model")) {
    stop("model must be returned by fit_crnp", call. = FALSE)
  }
  target_matrix <- .crnp_validate_matrix(target_matrix, "target_matrix")
  if (!all(model$observed_proteins %in% colnames(target_matrix))) {
    stop("target_matrix does not contain every measured model protein",
         call. = FALSE)
  }
  x <- target_matrix[, model$observed_proteins, drop = FALSE]
  standardized <- .crnp_standardize(x, center, scale)
  score <- model$calibration_intercept + model$calibration_slope *
    as.numeric(standardized$matrix %*% model$effective_weights)
  if (standardize_score) score <- .crnp_z_vector(score)
  attr(score, "target_center") <- standardized$center
  attr(score, "target_scale") <- standardized$scale
  score
}

#' Print a CRNP model
#' @param x A fitted CRNP model.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.crnp_model <- function(x, ...) {
  cat("Coverage-robust network projection model (v", x$version, ")\n", sep = "")
  cat("  Module proteins:", length(x$weights), "\n")
  cat("  Measured proteins:", length(x$observed_proteins), "\n")
  cat("  Structural coverage:", sprintf("%.1f%%", 100 * x$coverage_fraction), "\n")
  cat("  Recovery fraction:", sprintf("%.3f", x$recovery_fraction), "\n")
  cat("  Selected lambda:", format(x$lambda), "\n")
  cat("  Source-implied correlation fidelity:",
      sprintf("%.3f", x$correlation_fidelity), "\n")
  invisible(x)
}

