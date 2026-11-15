.crnp_z_vector <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

.crnp_validate_matrix <- function(x, name, require_finite = TRUE) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (!nrow(x) || !ncol(x)) stop(name, " must have at least one row and column", call. = FALSE)
  if (is.null(colnames(x)) || anyNA(colnames(x)) || anyDuplicated(colnames(x))) {
    stop(name, " must have unique, nonmissing column names", call. = FALSE)
  }
  if (require_finite && any(!is.finite(x))) {
    stop(name, " must contain only finite values", call. = FALSE)
  }
  x
}

.crnp_standardize <- function(x, center = NULL, scale = NULL) {
  x <- as.matrix(x)
  if (is.null(center)) center <- colMeans(x)
  if (is.null(scale)) scale <- apply(x, 2L, stats::sd)
  if (length(center) != ncol(x) || length(scale) != ncol(x)) {
    stop("center and scale must match the number of matrix columns", call. = FALSE)
  }
  scale[!is.finite(scale) | scale == 0] <- 1
  z <- sweep(sweep(x, 2L, center, "-"), 2L, scale, "/")
  list(matrix = z, center = center, scale = scale)
}

.crnp_robust_cor <- function(x) {
  r <- WGCNA::bicor(
    x, use = "pairwise.complete.obs", maxPOutliers = 0.05, quick = 0
  )
  r[!is.finite(r)] <- 0
  r <- (r + t(r)) / 2
  diag(r) <- 1
  r
}

.crnp_safe_solve <- function(a, b) {
  ans <- tryCatch(solve(a, b), error = function(e) NULL)
  if (is.null(ans)) ans <- qr.solve(a + diag(1e-8, nrow(a)), b, tol = 1e-10)
  as.numeric(ans)
}

.crnp_effective_weights <- function(r, w, observed, lambda, recovery_fraction = 1) {
  missing <- setdiff(seq_along(w), observed)
  if (!length(missing)) return(w[observed])
  roo <- r[observed, observed, drop = FALSE]
  rom <- r[observed, missing, drop = FALSE]
  correction <- .crnp_safe_solve(
    roo + diag(lambda, length(observed)), rom %*% w[missing]
  )
  w[observed] + recovery_fraction * correction
}

.crnp_recovery_fraction <- function(coverage_fraction) {
  max(0, min(1, coverage_fraction * (coverage_fraction - 0.20)))
}

.crnp_ridge_weights <- function(r, w, observed, lambda) {
  roo <- r[observed, observed, drop = FALSE]
  cov_oy <- as.numeric(r[observed, , drop = FALSE] %*% w)
  .crnp_safe_solve(roo + diag(lambda, length(observed)), cov_oy)
}

.crnp_pc1_weights <- function(r, w, observed) {
  roo <- r[observed, observed, drop = FALSE]
  ev <- eigen(roo, symmetric = TRUE, only.values = FALSE)$vectors[, 1L]
  if (sum(ev * w[observed]) < 0) ev <- -ev
  as.numeric(ev)
}

.crnp_topology_order <- function(r, power = 8) {
  adjacency <- ((1 + r) / 2)^power
  diag(adjacency) <- 1
  stats::hclust(stats::as.dist(1 - adjacency), method = "average")$order
}

.crnp_make_mask <- function(w, r, fraction, mechanism, seed, power = 8) {
  set.seed(seed)
  p <- length(w)
  k <- max(2L, min(p - 1L, as.integer(round(p * fraction))))
  aw <- abs(w)
  eps <- max(aw, na.rm = TRUE) * 0.02 + 1e-8
  if (mechanism == "random") {
    observed <- sample.int(p, k)
  } else if (mechanism == "hub_depleted") {
    observed <- sample.int(p, k, prob = (max(aw) - aw) + eps)
  } else if (mechanism == "hub_enriched") {
    observed <- sample.int(p, k, prob = aw + eps)
  } else if (mechanism == "topology_blocked") {
    ord <- .crnp_topology_order(r, power = power)
    n_missing <- p - k
    start <- sample.int(p, 1L)
    idx <- ((start - 1L + seq_len(n_missing) - 1L) %% p) + 1L
    observed <- setdiff(seq_len(p), ord[idx])
  } else {
    stop("Unknown coverage mechanism: ", mechanism, call. = FALSE)
  }
  sort(unique(observed))
}

.crnp_tune_regularization <- function(x, w, masks, lambdas,
                                      method = c("crnp", "ridge"), folds = 3L,
                                      recovery_fraction = 1) {
  method <- match.arg(method)
  set.seed(93017L + ncol(x))
  fold_id <- sample(rep(seq_len(folds), length.out = nrow(x)))
  y <- .crnp_z_vector(as.numeric(x %*% w))
  fold_r <- lapply(seq_len(folds), function(f) {
    .crnp_robust_cor(x[fold_id != f, , drop = FALSE])
  })
  detail <- vector("list", length(lambdas) * length(masks))
  pos <- 0L
  for (lambda in lambdas) {
    for (m in seq_along(masks)) {
      observed <- masks[[m]]
      pred <- rep(NA_real_, nrow(x))
      for (f in seq_len(folds)) {
        beta <- if (method == "crnp") {
          .crnp_effective_weights(
            fold_r[[f]], w, observed, lambda, recovery_fraction
          )
        } else {
          .crnp_ridge_weights(fold_r[[f]], w, observed, lambda)
        }
        pred[fold_id == f] <- as.numeric(
          x[fold_id == f, observed, drop = FALSE] %*% beta
        )
      }
      pos <- pos + 1L
      detail[[pos]] <- data.frame(
        lambda = lambda,
        mask = m,
        pearson_r = stats::cor(y, pred),
        standardized_mse = mean((y - .crnp_z_vector(pred))^2),
        stringsAsFactors = FALSE
      )
    }
  }
  detail <- do.call(rbind, detail)
  split_detail <- split(detail, detail$lambda)
  summary <- do.call(rbind, lapply(split_detail, function(d) {
    data.frame(
      lambda = d$lambda[[1L]],
      median_r = stats::median(d$pearson_r, na.rm = TRUE),
      lower_quartile_r = unname(stats::quantile(d$pearson_r, 0.25, na.rm = TRUE)),
      mean_mse = mean(d$standardized_mse, na.rm = TRUE)
    )
  }))
  summary <- summary[order(
    -summary$median_r, -summary$lower_quartile_r, summary$mean_mse,
    summary$lambda
  ), , drop = FALSE]
  rownames(summary) <- NULL
  list(best = summary$lambda[[1L]], summary = summary, detail = detail)
}

.crnp_calibration <- function(source_z, truth, beta, observed) {
  raw <- as.numeric(source_z[, observed, drop = FALSE] %*% beta)
  fit <- stats::lm(truth ~ raw)
  c(intercept = unname(stats::coef(fit)[[1L]]),
    slope = unname(stats::coef(fit)[[2L]]))
}

.crnp_method_weights <- function(method, r, w, observed, lambda_crnp,
                                 lambda_ridge, recovery_fraction) {
  if (method == "direct_kME") {
    factor <- sum(abs(w)) / sum(abs(w[observed]))
    return(w[observed] * factor)
  }
  if (method == "source_mean_imputed_fixed") return(w[observed])
  if (method == "signed_mean") return(sign(w[observed]) / length(observed))
  if (method == "subset_PC1") return(.crnp_pc1_weights(r, w, observed))
  if (method == "ridge") return(.crnp_ridge_weights(r, w, observed, lambda_ridge))
  if (method == "CRNP_full") {
    return(.crnp_effective_weights(r, w, observed, lambda_crnp, 1))
  }
  if (method == "CRNP") {
    return(.crnp_effective_weights(
      r, w, observed, lambda_crnp, recovery_fraction
    ))
  }
  stop("Unknown projection method: ", method, call. = FALSE)
}

.crnp_ccc <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3L) return(NA_real_)
  2 * stats::cov(x, y) /
    (stats::var(x) + stats::var(y) + (mean(x) - mean(y))^2)
}

.crnp_score_metrics <- function(reference, estimate) {
  keep <- is.finite(reference) & is.finite(estimate)
  reference <- reference[keep]
  estimate <- estimate[keep]
  if (length(reference) < 5L || stats::sd(estimate) == 0) {
    return(c(pearson_r = NA_real_, spearman_rho = NA_real_, ccc = NA_real_,
             standardized_rmse = NA_real_, calibration_intercept = NA_real_,
             calibration_slope = NA_real_, extreme_tertile_retention = NA_real_,
             normalized_rank_displacement = NA_real_))
  }
  fit <- stats::lm(reference ~ estimate)
  n <- length(reference)
  cut_low <- stats::quantile(reference, 1 / 3, names = FALSE)
  cut_high <- stats::quantile(reference, 2 / 3, names = FALSE)
  ref_group <- ifelse(reference <= cut_low, -1L,
                      ifelse(reference >= cut_high, 1L, 0L))
  est_low <- stats::quantile(estimate, 1 / 3, names = FALSE)
  est_high <- stats::quantile(estimate, 2 / 3, names = FALSE)
  est_group <- ifelse(estimate <= est_low, -1L,
                      ifelse(estimate >= est_high, 1L, 0L))
  extreme <- ref_group != 0L
  c(
    pearson_r = stats::cor(reference, estimate, method = "pearson"),
    spearman_rho = stats::cor(reference, estimate, method = "spearman"),
    ccc = .crnp_ccc(reference, estimate),
    standardized_rmse = sqrt(mean((.crnp_z_vector(reference) -
      .crnp_z_vector(estimate))^2)),
    calibration_intercept = unname(stats::coef(fit)[[1L]]),
    calibration_slope = unname(stats::coef(fit)[[2L]]),
    extreme_tertile_retention = mean(ref_group[extreme] == est_group[extreme]),
    normalized_rank_displacement = mean(abs(rank(reference) - rank(estimate))) /
      max(1, n - 1)
  )
}

