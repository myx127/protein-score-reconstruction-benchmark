#' Calculate source-implied fidelity of a partial projection
#'
#' Calculates the squared-correlation and variance-normalized error fidelity
#' implied by the source correlation matrix for the implemented coefficients.
#'
#' @param source_correlation Source protein correlation matrix.
#' @param weights Frozen signed weights, in the order of the correlation matrix.
#' @param observed_proteins Integer positions or protein names that are measured.
#' @param lambda Positive regularization parameter.
#' @param recovery_fraction Fraction of the source-implied correction retained.
#'
#' @return Named numeric vector containing `correlation_fidelity` and
#'   `mse_fidelity`.
#' @export
score_fidelity <- function(source_correlation, weights, observed_proteins,
                           lambda, recovery_fraction = 1) {
  r <- as.matrix(source_correlation)
  w <- as.numeric(weights)
  if (!is.numeric(observed_proteins)) {
    protein_names <- names(weights)
    if (is.null(protein_names)) protein_names <- colnames(r)
    observed <- match(observed_proteins, protein_names)
    if (anyNA(observed)) stop("Unknown observed protein", call. = FALSE)
  } else {
    observed <- as.integer(observed_proteins)
  }
  if (nrow(r) != length(w) || ncol(r) != length(w)) {
    stop("source_correlation and weights have incompatible dimensions",
         call. = FALSE)
  }
  b <- .crnp_effective_weights(r, w, observed, lambda, recovery_fraction)
  total_var <- as.numeric(crossprod(w, r %*% w))
  pred_var <- as.numeric(crossprod(
    b, r[observed, observed, drop = FALSE] %*% b
  ))
  cov_full_pred <- as.numeric(crossprod(
    b, r[observed, , drop = FALSE] %*% w
  ))
  if (!is.finite(total_var) || !is.finite(pred_var) ||
      total_var <= 0 || pred_var <= 0) {
    return(c(correlation_fidelity = NA_real_, mse_fidelity = NA_real_))
  }
  correlation_fidelity <- cov_full_pred^2 / (total_var * pred_var)
  mse <- total_var - 2 * cov_full_pred + pred_var
  c(
    correlation_fidelity = max(0, min(1, correlation_fidelity)),
    mse_fidelity = 1 - mse / total_var
  )
}

