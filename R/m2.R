#-------------------------------------------------------------------------------
# Compute M2 and M2* for cogirt IRT models
#
# Features:
# - NA-safe (pairwise complete)
# - Works for multidimensional, SDT, contrast-coded models
# - Diagonal-weight M2 (fast, mirt-style)
# - Robust M2* using df-adjusted correction (valid for diag-W)
#
# @export m2
#-------------------------------------------------------------------------------

m2 <- function(obj, use_pairs = TRUE) {

  #------------------------------------------------------------
  # 0. Extract data
  #------------------------------------------------------------
  Y <- obj$y
  Y <- as.matrix(Y)

  N <- nrow(Y)
  I <- ncol(Y)

  #------------------------------------------------------------
  # 1. Model-implied probabilities
  #------------------------------------------------------------
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

  P <- as.matrix(P)

  #------------------------------------------------------------
  # 2. Univariate observed margins (NA-safe)
  #------------------------------------------------------------
  N1 <- colSums(!is.na(Y))  # available N for each item

  p1_obs <- numeric(I)
  for (j in seq_len(I)) {
    if (N1[j] > 0L) {
      p1_obs[j] <- mean(Y[!is.na(Y[, j]), j])
    } else {
      p1_obs[j] <- NA_real_
    }
  }

  p1_mod <- colMeans(P)

  var1 <- p1_mod * (1 - p1_mod) / pmax(N1, 1L)
  var1[N1 <= 0L] <- NA_real_

  #------------------------------------------------------------
  # 3. Bivariate margins (pairwise complete)
  #------------------------------------------------------------
  p2_obs <- NULL
  p2_mod <- NULL
  var2   <- NULL
  N2     <- NULL

  if (use_pairs && I > 1L) {

    p2_obs_vec <- c()
    p2_mod_vec <- c()
    N2_vec     <- c()

    for (i in seq_len(I - 1L)) {
      for (j in seq.int(i + 1L, I)) {

        cc <- !is.na(Y[, i]) & !is.na(Y[, j])
        Nij <- sum(cc)

        if (Nij > 0L) {
          p2_obs_ij <- mean(Y[cc, i] * Y[cc, j])
          p2_mod_ij <- mean(P[cc, i] * P[cc, j])

          p2_obs_vec <- c(p2_obs_vec, p2_obs_ij)
          p2_mod_vec <- c(p2_mod_vec, p2_mod_ij)
          N2_vec     <- c(N2_vec, Nij)
        } else {
          p2_obs_vec <- c(p2_obs_vec, NA_real_)
          p2_mod_vec <- c(p2_mod_vec, NA_real_)
          N2_vec     <- c(N2_vec, 0L)
        }
      }
    }

    p2_obs <- p2_obs_vec
    p2_mod <- p2_mod_vec
    N2     <- N2_vec

    var2 <- p2_mod * (1 - p2_mod) / pmax(N2, 1L)
    var2[N2 <= 0L] <- NA_real_
  }

  #------------------------------------------------------------
  # 4. Stack margins
  #------------------------------------------------------------
  s_obs <- p1_obs
  s_mod <- p1_mod
  var_all <- var1

  if (!is.null(p2_obs)) {
    s_obs   <- c(s_obs, p2_obs)
    s_mod   <- c(s_mod, p2_mod)
    var_all <- c(var_all, var2)
  }

  d <- s_obs - s_mod

  #------------------------------------------------------------
  # 5. Remove unusable margins
  #------------------------------------------------------------
  eps <- 1e-10
  bad <- is.na(var_all) | (var_all < eps)
  keep <- !bad

  if (all(bad)) {
    stop("All margins have zero or undefined variance; M2 cannot be computed.")
  }

  d_use   <- d[keep]
  var_use <- var_all[keep]

  #------------------------------------------------------------
  # 6. Compute diagonal-weight M2
  #------------------------------------------------------------
  W_inv <- 1 / var_use
  M2 <- sum(d_use^2 * W_inv)

  #------------------------------------------------------------
  # 7. Degrees of freedom
  #------------------------------------------------------------
  n_margins_kept <- sum(keep)
  n_par <- if (!is.null(obj$par)) obj$par else I

  df <- n_margins_kept - n_par
  if (df < 1L) df <- 1L

  p_val <- stats::pchisq(M2, df = df, lower.tail = FALSE)

  #------------------------------------------------------------
  # 8. RMSEA and SRMSR
  #------------------------------------------------------------
  num <- M2 - df
  den <- df * (N - 1)

  RMSEA_M2 <- sqrt(max(num / den, 0))
  SRMSR_M2 <- sqrt(mean(d_use^2))

  #------------------------------------------------------------
  # 9. Robust M2* using df-adjustment (Cai & Hansen 2013, diag-W form)
  #------------------------------------------------------------
  k_scale <- n_margins_kept / df  # scaling > 1

  df_star <- df * k_scale
  M2_star <- M2
  p_star  <- stats::pchisq(M2_star, df = df_star, lower.tail = FALSE)

  #------------------------------------------------------------
  # 10. Return results
  #------------------------------------------------------------
  list(
    type     = "Fast M2 (diag W, NA-safe)",
    M2       = M2,
    df       = df,
    p        = p_val,
    RMSEA_M2 = RMSEA_M2,
    SRMSR_M2 = SRMSR_M2,

    M2_star  = M2_star,
    df_star  = df_star,
    p_star   = p_star,

    n_margins = n_margins_kept,
    n_par     = n_par,

    details = list(
      link      = obj$link,
      use_pairs = use_pairs,
      N         = N,
      I         = I
    )
  )
}
