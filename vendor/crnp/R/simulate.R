#' Simulate a fixed protein module for examples and tests
#'
#' Generates phenotype-free source and target protein matrices from a block
#' correlation model. The function is provided for reproducible examples; it
#' is not part of CRNP estimation.
#'
#' @param n_source Number of source participants.
#' @param n_target Number of target participants.
#' @param proteins Number of module proteins.
#' @param blocks Number of correlated protein blocks.
#' @param shift Strength of target covariance shift.
#' @param seed Random seed.
#'
#' @return List with source matrix, target matrix, and named weights.
#' @export
simulate_crnp_data <- function(n_source = 160L, n_target = 35L,
                               proteins = 30L, blocks = 3L,
                               shift = 0.10, seed = 1L) {
  set.seed(seed)
  proteins <- as.integer(proteins)
  block <- rep(seq_len(blocks), length.out = proteins)
  sigma <- matrix(0.08, proteins, proteins)
  sigma[outer(block, block, "==")] <- 0.55
  diag(sigma) <- 1
  sigma_target <- (1 - shift) * sigma + shift * diag(proteins)
  draw <- function(n, s) {
    eig <- eigen(s, symmetric = TRUE)
    root <- eig$vectors %*% diag(sqrt(pmax(eig$values, 1e-8)), proteins)
    matrix(stats::rnorm(n * proteins), n, proteins) %*% t(root)
  }
  protein_names <- sprintf("P%03d", seq_len(proteins))
  source <- draw(n_source, sigma)
  target <- draw(n_target, sigma_target)
  colnames(source) <- colnames(target) <- protein_names
  rownames(source) <- sprintf("S%04d", seq_len(n_source))
  rownames(target) <- sprintf("T%03d", seq_len(n_target))
  weights <- stats::setNames(seq(1.2, 0.25, length.out = proteins), protein_names)
  weights[seq(4L, proteins, by = 5L)] <- -weights[seq(4L, proteins, by = 5L)]
  list(source_matrix = source, target_matrix = target, weights = weights)
}

