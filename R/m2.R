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

m2 <- function(object,
                            n_mc   = 5000L,
                            link   = c("logit","probit"),
                            weight = c("diag","full"),
                            n_par  = NULL,
                            eps    = 1e-10) {

  link   <- match.arg(link)
  weight <- match.arg(weight)

  Y <- object$y
  if (is.null(Y)) stop("object$y is missing.", call. = FALSE)
  N <- nrow(Y); I <- ncol(Y)

  ## ----- 1. population latent draws (NOT person scores) -----
  mu  <- object$omega_mu
  Sig <- object$omega_sigma2
  if (is.null(mu)) {
    d   <- 1L
    mu  <- 0
    Sig <- 1
  } else {
    d <- length(mu)
  }

  if (d == 1L) {
    TH <- matrix(stats::rnorm(n_mc, mu, sqrt(Sig)), ncol = 1L)
  } else {
    if (!requireNamespace("mvtnorm", quietly = TRUE))
      stop("mvtnorm required.", call. = FALSE)
    TH <- mvtnorm::rmvnorm(n_mc, mean = mu, sigma = Sig)
  }
  wts <- rep(1 / n_mc, n_mc)

  ## ----- 2. item-level params only -----
  # item intercepts
  nu <- object$nu1
  if (is.null(nu)) nu <- rep(0, I)

  # item slopes
  lambda <- object$lambda1
  if (is.null(lambda)) {
    lambda <- matrix(1, nrow = I, ncol = d)
  } else {
    # make it I x d
    if (is.null(dim(lambda))) {
      lambda <- matrix(lambda, nrow = I, ncol = 1L)
    } else if (nrow(lambda) == d && ncol(lambda) == I) {
      lambda <- t(lambda)
    }
  }

  ## IMPORTANT: drop all person-level / condition-level stuff
  # gamma0, zeta0, kappa0 are NOT used here

  ## ----- 3. model-implied probs from population model -----
  # TH: n_mc x d
  # lambda: I x d
  eta <- TH %*% t(lambda)                # n_mc x I
  eta <- sweep(eta, 2, nu, "+")          # add intercepts

  if (link == "logit") {
    P <- 1 / (1 + exp(-eta))
  } else {
    P <- pnorm(eta)
  }

  ## ----- 4. model-implied margins -----
  # 1-way
  p1 <- as.numeric(colSums(P * wts))

  # 2-way
  p2 <- matrix(0, I, I)
  for (i in 1:I) {
    pi <- P[, i]
    for (j in i:I) {
      val <- sum(wts * pi * P[, j])
      p2[i, j] <- val
      p2[j, i] <- val
    }
  }

  # 3-way
  p3 <- array(0, dim = c(I, I, I))
  for (i in 1:I) {
    Pi <- P[, i]
    for (j in i:I) {
      Pj <- P[, j]
      for (k in j:I) {
        val <- sum(wts * Pi * Pj * P[, k])
        p3[i, j, k] <- val; p3[i, k, j] <- val
        p3[j, i, k] <- val; p3[j, k, i] <- val
        p3[k, i, j] <- val; p3[k, j, i] <- val
      }
    }
  }

  # 4-way on demand
  p4_env <- new.env(parent = emptyenv())
  get_p4 <- function(a,b,c,d) {
    key <- paste(sort(c(a,b,c,d)), collapse = ".")
    if (!exists(key, envir = p4_env, inherits = FALSE)) {
      val <- sum(wts * P[,a] * P[,b] * P[,c] * P[,d])
      assign(key, val, envir = p4_env)
    }
    get(key, envir = p4_env, inherits = FALSE)
  }

  ## ----- 5. observed margins -----
  obs1 <- colMeans(Y)
  obs2 <- crossprod(Y) / N

  pairs   <- which(upper.tri(matrix(0, I, I)), arr.ind = TRUE)
  s_obs   <- c(obs1, obs2[upper.tri(obs2, diag = FALSE)])
  s_mod   <- c(p1,   p2[upper.tri(p2,   diag = FALSE)])
  diff    <- s_obs - s_mod
  p_stats <- length(diff)

  ## ----- 6. Xi (Cai duplicate-aware), but from CLEAN margins -----
  Xi <- matrix(0, p_stats, p_stats)

  # uni-uni
  for (i in 1:I) {
    for (j in i:I) {
      cij <- p2[i, j] - p1[i] * p1[j]
      Xi[i, j] <- cij
      Xi[j, i] <- cij
    }
  }

  # uni-bi
  for (r in 1:I) {
    for (ab in seq_len(nrow(pairs))) {
      i <- pairs[ab, 1]; j <- pairs[ab, 2]
      crij <- p3[r, i, j] -
        p1[r] * p2[i, j] -
        p1[i] * p2[r, j] -
        p1[j] * p2[r, i] +
        2 * p1[r] * p1[i] * p1[j]
      Xi[r, I + ab] <- crij
      Xi[I + ab, r] <- crij
    }
  }

  # bi-bi
  for (ab in seq_len(nrow(pairs))) {
    i <- pairs[ab, 1]; j <- pairs[ab, 2]
    for (cd in ab: nrow(pairs)) {
      k <- pairs[cd, 1]; l <- pairs[cd, 2]

      if (i == k && j == l) {
        cov_val <- p2[i, j] - p2[i, j]^2
      } else {
        idxs <- c(i, j, k, l)
        u <- length(unique(idxs))
        if (u == 4L) {
          pijkl  <- get_p4(i, j, k, l)
          cov_val <- pijkl - p2[i, j] * p2[k, l]
        } else if (u == 3L) {
          common <- intersect(c(i, j), c(k, l))[1]
          others <- setdiff(idxs, common)
          o1 <- others[1]; o2 <- others[2]
          cov_val <- p3[common, o1, o2] - p2[common, o1] * p2[common, o2]
        } else {
          cov_val <- p2[i, j] - p2[i, j]^2
        }
      }

      Xi[I + ab, I + cd] <- cov_val
      Xi[I + cd, I + ab] <- cov_val
    }
  }

  Sigma <- Xi / N

  ## ----- 7. weight -----
  if (weight == "diag") {
    W <- diag(1 / (diag(Sigma) + eps), nrow = p_stats)
  } else {
    W <- tryCatch(
      solve(Sigma + eps * diag(p_stats)),
      error = function(e) diag(1 / (diag(Sigma) + eps), nrow = p_stats)
    )
  }

  ## ----- 8. parameter count (mirt-like) -----
  if (is.null(n_par)) {
    # special case: your fit is 1p / Rasch-like
    mdl <- if (!is.null(object$model)) tolower(object$model) else ""
    if (mdl %in% c("1p","rasch","1pl")) {
      n_par <- I + 1L      # mirt: 10 items -> 11 pars -> 55-11=44 df
    } else {
      # generic: item intercepts + item slopes + latent mean + latent cov
      n_par <- length(nu) + length(lambda) + d + d*(d+1)/2
    }
  }

  df <- p_stats - n_par
  if (df < 1L) df <- 1L

  ## ----- 9. test -----
  M2    <- as.numeric(N * crossprod(diff, W %*% diff))
  pval  <- stats::pchisq(M2, df = df, lower.tail = FALSE)
  RMSEA <- sqrt(max((M2 - df) / (N * df), 0))

  list(
    type  = "M2",
    stat  = M2,
    df    = df,
    p     = pval,
    RMSEA = RMSEA,
    details = list(
      integ   = "mc",
      n_mc    = n_mc,
      weight  = weight,
      n_stats = p_stats,
      n_par   = n_par
    )
  )
}
