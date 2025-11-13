#-------------------------------------------------------------------------------
#' Yen's Q3
#'
#' Computes Yen's Q3 (residual correlations) and associated posterior predictive
#' p-values (PPP) computed via parametric replication.
#
#' @param object An object of class 'cog_irt'.
#' @param num_ppc Number of posterior predictive replications. Default 200.
#
#' @return A list with elements for fit information q3, q3_avg, q3_max,
#' q3_avg_ppp, q3_max_ppp, q3_pair_ppp
#
#' @references:
#' Yen, W. M. (1984). Effects of local item dependence on the fit and equating
#' performance of the three-parameter logistic model. \emph{Applied
#' Psychological Measurement, 8}, 125-145.
#'
#' @export q3
#-------------------------------------------------------------------------------

q3 <- function(object, num_ppc = 200) {

  obj <- object
  Y <- obj$y
  N <- nrow(Y)
  I <- ncol(Y)

  # Model-implied probabilities (P_ij)
  P <- dich_response_model(y = obj$y, omega = obj$omega1, gamma = obj$gamma0,
                           lambda = obj$lambda1, zeta = obj$zeta0, nu = obj$nu1,
                           kappa = obj$kappa0, link = obj$link)$p

  # Q3 Statistics --------------------------------------------------------------
  q3_func <- function(Y, P){
    U  <- Y - P
    rm <- rowMeans(U, na.rm = TRUE)      # person means
    Uc <- U - rm                         # person-centered residuals
    v  <- apply(Uc, 2, var, na.rm = TRUE)
    keep <- is.finite(v) & v > 0
    R <- matrix(NA_real_, ncol(Uc), ncol(Uc))
    if (any(keep)) {
      Rk <- stats::cor(Uc[, keep, drop = FALSE], use = "pairwise.complete.obs")
      R[keep, keep] <- Rk
    }
    diag(R) <- NA_real_
    R
  }

  q3 <- q3_func(Y, P)

  q3_avg <- mean(q3, na.rm = TRUE)
  q3_max <- max(q3, na.rm = TRUE)

  # Posterior predictive p-values ----------------------------------------------
  q3_avg_rep <- numeric(0)
  q3_max_rep <- numeric(0)
  q3_pair_ppp <- NULL

  for (b in seq_len(num_ppc)) {
    Yrep <- matrix(stats::rbinom(N * I, 1, prob = as.vector(P)), nrow = N)
    q3r <- q3_func(Yrep, P)
    q3_avg_rep <- c(q3_avg_rep, mean(q3r, na.rm = TRUE))
    q3_max_rep <- c(q3_max_rep, max(q3r,  na.rm = TRUE))
  }

  q3_avg_ppp <- mean(q3_avg_rep >= q3_avg, na.rm = TRUE)
  q3_max_ppp <- mean(q3_max_rep >= q3_max, na.rm = TRUE)

  lowtri_idx <- which(lower.tri(q3))
  rep_mat <- matrix(NA_real_, num_ppc, length(lowtri_idx))
  for (b in seq_len(num_ppc)) {
    Yrep <- matrix(stats::rbinom(N * I, 1, as.vector(P)), nrow = N)
    q3r <- q3_func(Yrep, P)
    rep_mat[b, ] <- q3r[lowtri_idx]
  }
  vobs <- q3[lowtri_idx]
  pair_ppp <- colMeans(t(t(rep_mat) >= vobs), na.rm = TRUE)
  q3_pair_ppp <- matrix(NA_real_, I, I)
  diag(q3_pair_ppp) <- NA_real_
  q3_pair_ppp[lowtri_idx] <- pair_ppp

  list(
    q3 = q3,
    q3_avg = q3_avg,
    q3_max = q3_max,
    q3_avg_ppp = q3_avg_ppp,
    q3_max_ppp = q3_max_ppp,
    q3_pair_ppp = q3_pair_ppp
  )
}
