#-------------------------------------------------------------------------------
#' Compute M2 and M2* for IRT model object via Monte Carlo
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
#' Cai, L. & Hansen, M. (2013). Limited-information goodness-of-object testing of
#' hierarchical item factor models. \emph{British Journal of Mathematical and
#' Statistical Psychology, 66}, 245-276.
#'
#' @export m2
#-------------------------------------------------------------------------------

m2 <- function(
    obj,
    use_pairs = TRUE
) {
  Y <- obj$y
  N <- nrow(Y)
  I <- ncol(Y)

  # model-implied probs, person-by-item
  P <- dich_response_model(
    y      = obj$y,
    omega  = obj$omega1,
    gamma  = obj$gamma0,
    lambda = obj$lambda1,
    zeta   = obj$zeta0,
    nu     = obj$nu1,
    kappa  = obj$kappa0,
    link   = obj$link
  )$p

  # 1-way observed and model
  p1_obs <- colMeans(Y, na.rm = TRUE)
  p1_mod <- colMeans(P, na.rm = TRUE)

  # variances for 1-way (from 2^n table logic)
  var1 <- p1_mod * (1 - p1_mod) / N

  # stack
  s_obs <- p1_obs
  s_mod <- p1_mod
  w_diag <- 1 / pmax(var1, 1e-10)

  if (use_pairs && I > 1L) {
    # observed 2-way
    p2_obs <- (t(Y) %*% Y) / N
    # model 2-way: average of person-specific products
    p2_mod <- (t(P) %*% P) / N

    # we only want lower triangle, i < j
    idx <- which(lower.tri(p2_obs), arr.ind = TRUE)

    p2_obs_vec <- p2_obs[idx]
    p2_mod_vec <- p2_mod[idx]

    # var for 2-way margins (from 2^n tables): p(ij)*(1 - p(ij)) / N
    var2 <- p2_mod_vec * (1 - p2_mod_vec) / N

    s_obs <- c(s_obs, p2_obs_vec)
    s_mod <- c(s_mod, p2_mod_vec)
    w_diag <- c(w_diag, 1 / pmax(var2, 1e-10))
  }

  d <- s_obs - s_mod

  # diagonal-weight chi-square
  stat <- sum(d^2 * w_diag)

  # crude df: number of margins - number of free item params
  # you may want to plug in your own param counter here
  n_margins <- length(d)
  # placeholder: try to read model type
  model_name <- if (!is.null(obj$model)) obj$model else "2p"
  if (model_name %in% c("1p","rasch")) {
    n_par <- I  # one per item
  } else if (model_name == "2p") {
    n_par <- 2L * I
  } else if (model_name == "3p") {
    n_par <- 3L * I
  } else {
    n_par <- I
  }
  df <- max(n_margins - n_par, 1L)

  pval <- stats::pchisq(stat, df = df, lower.tail = FALSE)

  num <- stat - df
  den <- df * (N - 1)
  rmsea <- sqrt( max(num / den, 0) )


  list(
    type     = "MOJ-limited (diag)",
    stat     = stat,
    df       = df,
    p        = pval,
    rmsea = rmsea,
    n_marg   = n_margins,
    details  = list(
      link      = link,
      use_pairs = use_pairs,
      N         = N
    )
  )
}
