#-------------------------------------------------------------------------------
#' Compute M2 and M2* for IRT model fit via Monte Carlo
#'
#' Uses Monte Carlo draws from the latent trait distribution to estimate
#' model-implied 1- and 2-way margins, compares them to observed margins, and
#' returns M2, df, p-value, and an RMSEA-like index (M2*).
#'
#' @param object An object of class 'cog_irt'.
#' @param type "M2" or "M2*" (small-sample adjusted)
#' @param n_mc umber of Monte Carlo samples (draws) used for integration.
#' @return list with statistic, df, p, RMSEA2, and type
#'
#' @references
#'
#' Cai, L. & Hansen, M. (2013). Limited-information goodness-of-fit testing of
#' hierarchical item factor models. \emph{British Journal of Mathematical and
#' Statistical Psychology, 66}, 245-276.
#'
#' @export
#-------------------------------------------------------------------------------

m2 <- function(fit, type = "M2", link = c("probit","logit"), n_mc   = 50000,
               eps    = 1e-8) {
  if (is.null(fit$omega1)) {
    stop("Model does not appear to be fitted to data.",
         call. = FALSE)
  }
  K <- nrow(fit$y)
  IJ <- ncol(fit$y)

  OMEGA <- mvtnorm::rmvnorm(n_mc, mean = fit$omega_mu, sigma = fit$omega_sigma2)
  wts   <- rep(x = 1/n_mc, n_mc)
  Q <- nrow(x = OMEGA)

  # Compute model-implied probabilities at nodes
  P <- dich_response_model(
    y      = matrix(0, nrow(OMEGA), ncol(fit$y)),
    omega  = OMEGA,
    gamma  = fit$gamma0,
    lambda = fit$lambda1,
    zeta   = matrix(0, nrow(OMEGA), ncol(fit$zeta0)),
    nu     = fit$nu1,
    kappa  = fit$kappa0,
    link   = fit$link
  )$p

  # Model-implied 1- and 2-way margins
  p_item_hat  <- as.vector(crossprod(wts, P))
  pairs   <- utils::combn(IJ, 2L)
  p_pair_hat <- apply(pairs, 2L, function(ix) sum(P[, ix[1]] * P[, ix[2]] * wts))

  # Observed margins
  p_item_obs  <- colMeans(fit$y)
  p_pair_obs <- apply(pairs, 2L, function(ix) mean(fit$y[, ix[1]] * fit$y[, ix[2]]))

  # Discrepancies & working covariance
  d <- c(p_item_obs - p_item_hat, p_pair_obs - p_pair_hat)
  var_item  <- pmax(p_item_hat * (1 - p_item_hat) / K, eps)
  var_pair <- pmax(p_pair_hat * (1 - p_pair_hat) / K, eps)

  M2 <- sum((d^2) / c(var_item, var_pair))
  c1 <- (K - 1) / K
  c2 <- df / (K - 1)
  M2s <- c1 * (M2 - c2)

  p_count <- fit$par
  df <- max(length(d) - p_count, 1L)
  pval   <- stats::pchisq(ifelse(type == "M2", M2, M2s), df, lower.tail = FALSE)
  RMSEA2 <- sqrt(x = max((ifelse(type == "M2", M2, M2s) - df) / (df * (K - 1)), 0))
  SRMR_item = sqrt(mean((p_item_obs - p_item_hat)^2))
  SRMR_pair = sqrt(mean((p_pair_obs - p_pair_hat)^2))

  # ---- SRMSR (mirt-style): RMS of residual correlations ----
  build_corr_from_margins <- function(p1, p2, J, pairs, eps = 1e-12) {
    # p1: length-J item means           (E[Y_j])
    # p2: length-J(J-1)/2 joint means   (E[Y_j Y_k]) in same order as `pairs`
    # returns: J x J correlation matrix on the 0/1 scale

    # variances on the diagonal
    var_j <- pmax(p1 * (1 - p1), eps)
    Sigma <- matrix(0, J, J)
    diag(Sigma) <- var_j

    # fill covariances from joint probs
    # Cov(Y_j, Y_k) = E[Y_j Y_k] - E[Y_j] E[Y_k]
    for (col in seq_len(ncol(pairs))) {
      j <- pairs[1L, col]; k <- pairs[2L, col]
      cov_jk <- p2[col] - p1[j] * p1[k]
      Sigma[j, k] <- cov_jk
      Sigma[k, j] <- cov_jk
    }

    # convert to correlations (guard zero-variance)
    sd_j <- sqrt(var_j)
    denom <- outer(sd_j, sd_j, "*")
    R <- Sigma / denom
    diag(R) <- 1
    R[!is.finite(R)] <- NA_real_   # if any degenerate items
    R
  }

  # Observed and model-implied correlation matrices
  R_obs <- build_corr_from_margins(pj_obs,  pjk_obs,  J, pairs, eps)
  R_hat <- build_corr_from_margins(pj_hat,  pjk_hat,  J, pairs, eps)

  # SRMSR = RMS of off-diagonal residual correlations
  lt <- lower.tri(R_obs)
  resid_lower <- (R_obs - R_hat)[lt]
  SRMSR <- sqrt(mean(resid_lower^2, na.rm = TRUE))


  list(
    type = type,
    stat = ifelse(type == "M2", M2, M2s),
    df = df,
    p = pval,
    RMSEA2 = RMSEA2,
    SRMR_item = SRMR_item,
    SRMR_pair = SRMR_pair,
    SRMSR = SRMSR
  )
}



# fit1p <- cog_irt(
#   data = ex1$y, model = "1p", link = "probit"
# )
# m2(fit1p, type = "M2",  link = "probit")
# m2(fit1p, type = "M2*", link = "probit")
#
# fit2p <- cog_irt(
#   data = ex1$y, model = "2p", link = "probit"
# )
#
# m2(fit2p, type = "M2",  link = "probit")
# m2(fit2p, type = "M2*", link = "probit")
#
#
# fitsdt <- cog_irt(
#   data = ex2$y, model = "sdt", key = ex2$key, link = "probit"
# )
# m2(fitsdt, link = "probit")
# m2(fitsdt, type = "M2*", link = "probit")

#
# nback_fit_contr <- cog_irt(data = nback$y, model = "sdt",
                           # contrast_codes = "contr.poly", key = nback$key,
                           # num_conditions = length(unique(nback$condition)),
                           # num_contrasts = 2)
# m2(nback_fit_contr, link = "probit")
# m2(nback_fit_contr, type = "M2*", link = "probit")
