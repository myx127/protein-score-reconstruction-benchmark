#' Compare CRNP with frozen projection methods
#'
#' Applies CRNP and the six comparators frozen for external evaluation. Every
#' method uses source-only coefficients and affine calibration.
#'
#' @param model A fitted `crnp_model`.
#' @param target_matrix Numeric target matrix containing measured proteins.
#' @param methods Projection methods to return.
#' @param center Optional target reference means.
#' @param scale Optional target reference standard deviations.
#' @param standardize_scores Whether to standardize each method score.
#'
#' @return Long data frame with one row per participant and method.
#' @export
compare_methods <- function(model, target_matrix,
                            methods = c(
                              "direct_kME", "signed_mean", "subset_PC1",
                              "ridge", "CRNP_full", "CRNP",
                              "source_mean_imputed_fixed"
                            ),
                            center = NULL, scale = NULL,
                            standardize_scores = TRUE) {
  if (!inherits(model, "crnp_model")) {
    stop("model must be returned by fit_crnp", call. = FALSE)
  }
  unknown <- setdiff(methods, names(model$method_models))
  if (length(unknown)) stop("Unknown method: ", unknown[[1L]], call. = FALSE)
  target_matrix <- .crnp_validate_matrix(target_matrix, "target_matrix")
  if (!all(model$observed_proteins %in% colnames(target_matrix))) {
    stop("target_matrix lacks a structurally measured model protein", call. = FALSE)
  }
  x <- target_matrix[, model$observed_proteins, drop = FALSE]
  standardized <- .crnp_standardize(x, center, scale)
  ids <- rownames(x)
  if (is.null(ids)) ids <- sprintf("sample_%04d", seq_len(nrow(x)))
  results <- lapply(methods, function(method) {
    component <- model$method_models[[method]]
    score <- component$calibration[["intercept"]] +
      component$calibration[["slope"]] *
      as.numeric(standardized$matrix %*% component$effective_weights)
    if (standardize_scores) score <- .crnp_z_vector(score)
    data.frame(sample_id = ids, method = method, score = score,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, results)
}

#' Evaluate known-answer recovery under structural coverage loss
#'
#' Artificially masks proteins from a complete target matrix, reconstructs the
#' frozen complete score, and returns method-specific fidelity metrics.
#'
#' @param source_matrix Complete source matrix.
#' @param target_matrix Complete target matrix used to define the known answer.
#' @param weights Named frozen weights.
#' @param masks Named list of observed-protein character vectors.
#' @param methods Methods passed to [compare_methods()].
#' @param lambdas Candidate source-only regularization parameters.
#' @param seed Source mask seed.
#'
#' @return Data frame with one row per mask and method.
#' @export
evaluate_coverage <- function(source_matrix, target_matrix, weights, masks,
                              methods = c(
                                "direct_kME", "signed_mean", "subset_PC1",
                                "ridge", "CRNP_full", "CRNP",
                                "source_mean_imputed_fixed"
                              ),
                              lambdas = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10),
                              seed = 93017L) {
  source_matrix <- .crnp_validate_matrix(source_matrix, "source_matrix")
  target_matrix <- .crnp_validate_matrix(target_matrix, "target_matrix")
  if (!all(names(weights) %in% colnames(target_matrix))) {
    stop("target_matrix must contain the complete weighted module", call. = FALSE)
  }
  target_z <- .crnp_standardize(
    target_matrix[, names(weights), drop = FALSE]
  )$matrix
  reference <- .crnp_z_vector(as.numeric(target_z %*% weights))
  if (is.null(names(masks))) names(masks) <- paste0("mask_", seq_along(masks))
  results <- vector("list", length(masks))
  for (i in seq_along(masks)) {
    observed <- names(weights)[names(weights) %in% masks[[i]]]
    model <- fit_crnp(
      source_matrix, weights, observed, lambdas = lambdas,
      calibration_mechanisms = c("random", "hub_depleted", "hub_enriched",
                                 "topology_blocked"),
      calibration_replicates = 2L, folds = 3L, seed = seed
    )
    scores <- compare_methods(
      model, target_matrix[, observed, drop = FALSE], methods = methods,
      standardize_scores = TRUE
    )
    by_method <- split(scores, scores$method)
    rows <- lapply(names(by_method), function(method) {
      metrics <- .crnp_score_metrics(reference, by_method[[method]]$score)
      data.frame(
        mask_id = names(masks)[[i]],
        method = method,
        module_size = length(weights),
        observed_proteins = length(observed),
        coverage_fraction = length(observed) / length(weights),
        t(metrics), check.names = FALSE, stringsAsFactors = FALSE
      )
    })
    results[[i]] <- do.call(rbind, rows)
  }
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}

