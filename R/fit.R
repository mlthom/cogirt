#-------------------------------------------------------------------------------
#' Item fit
#'
#' Computes Yen's Q3 (residual correlations) and associated posterior predictive
#' p-values (PPP).
#
#' @param object An object of class 'cog_irt'.
#
#' @return A list with elements for fit information q3, ppp_q3
#
#' @references:
#' Yen, W. M. (1984). Effects of local item dependence on the fit and equating
#' performance of the three-parameter logistic model. \emph{Applied
#' Psychological Measurement, 8}, 125-145.
#'
#' Li, T., Xie, C., & Jiao, H. (2017). Assessing fit of alternative
#' unidimensional polytomous IRT models using posterior predictive model
#' checking. Psychological Methods, 22(2), 397-408.
#'
#' @export fit
#-------------------------------------------------------------------------------

fit <- function(object, ...) {

  ellipsis <- list(...)

  if (is.null(x = ellipsis$chains)) {
    chains <- 3
  } else {
    chains <- ellipsis$chains
  }
  if (is.null(x = ellipsis$burn)) {
    burn <- 100
  } else {
    burn <- ellipsis$burn
  }
  if (is.null(x = ellipsis$thin)) {
    thin <- 10
  } else {
    thin <- ellipsis$thin
  }
  if (is.null(x = ellipsis$min_tune)) {
    min_tune <- 25
  } else {
    min_tune <- ellipsis$min_tune
  }
  if (is.null(x = ellipsis$tune_int)) {
    tune_int <- 25
  } else {
    tune_int <- ellipsis$tune_int
  }
  if (is.null(x = ellipsis$max_tune)) {
    max_tune <- 100
  } else {
    max_tune <- ellipsis$max_tune
  }
  if (is.null(x = ellipsis$niter)) {
    niter <- 1100
  } else {
    niter <- ellipsis$niter
  }

  # Draw parameter sets from posterior distribution ----------------------------

  draws_per_chain <- ((niter - burn) / thin)
  total_draws <- chains * draws_per_chain

  mc_draws <- mhmc_mc(
    chains = chains, y = object$y, obj_fun = dich_response_model,
    link = object$link, est_omega = object$est_omega,
    est_lambda = object$est_lambda, est_zeta = object$est_zeta,
    est_nu = object$est_nu, omega0 = object$omega1, gamma0 = object$gamma0,
    lambda0 = object$lambda1, zeta0 = object$zeta0, nu0 = object$nu1,
    kappa0 = object$kappa0, omega_mu = object$omega_mu,
    omega_sigma2 = object$omega_sigma2, lambda_mu = object$lambda_mu,
    lambda_sigma2 = object$lambda_sigma2, zeta_mu = object$zeta_mu,
    zeta_sigma2 = object$zeta_sigma2, nu_mu = object$nu_mu,
    nu_sigma2 = object$nu_sigma2, burn = burn, thin = thin, min_tune = min_tune,
    tune_int = tune_int, max_tune = max_tune, niter = niter, psrf = TRUE
  )

  # Created replicated dataset predicted values  -------------------------------

  rep_p <- vector("list", length = total_draws)
  for (chain in 1:chains) {
    for(draw in 1:((niter - burn) / thin)) {
      rep_p[[ (chain - 1) * draws_per_chain + draw]] <- dich_response_model(
        y = object$y,
        omega = if (object$est_omega) {
          matrix(
            data = mc_draws$mhmcdraws[[chain]]$omega_draws[draw, , ],
            nrow = nrow(object$omega1),
            ncol = ncol(object$omega1),
            byrow = TRUE
          )
        } else {
          object$omega1
        },
        gamma = object$gamma0,
        lambda = if (object$est_lambda) {
          matrix(
            data = mc_draws$mhmcdraws[[chain]]$lambda_draws[draw, , ],
            nrow = nrow(object$lambda1),
            ncol = ncol(object$lambda1),
            byrow = TRUE
          )
        } else {
          object$lambda1
        },
        zeta = if (object$est_zeta) {
          matrix(
            data = mc_draws$mhmcdraws[[chain]]$zeta_draws[draw, , ],
            nrow = nrow(object$zeta0),
            ncol = ncol(object$zeta0),
            byrow = TRUE
          )
        } else {
          object$zeta0
        },
        nu = if (object$est_nu) {
          matrix(
            data = mc_draws$mhmcdraws[[chain]]$nu_draws[draw, , ],
            nrow = nrow(object$nu1),
            ncol = ncol(object$nu1),
            byrow = TRUE
          )
        } else {
          object$nu1
        },
        kappa = object$kappa0,
        link = object$link
      )$p
    }
  }

  # Created replicated dataset observed values  --------------------------------

  rep_y <- lapply(
    X = rep_p,
    FUN = function(x) {
      y_vec <- rbinom(n = length(x), size = 1, prob = as.vector(x))
      matrix(y_vec, nrow = nrow(x), ncol = ncol(x))
    }
  )


  # Compute predicted values for observed data ---------------------------------

  p <- dich_response_model(
    y = object$y, omega = object$omega1, gamma = object$gamma0,
    lambda = object$lambda1, zeta = object$zeta0, nu = object$nu1,
    kappa = object$kappa0, link = object$link
  )$p

  # Q3 function ----------------------------------------------------------------
  q3_func <- function(y, p){
    r  <- y - p
    rmat <- stats::cor(r, use = "pairwise.complete.obs")
    diag(rmat) <- NA_real_
    rmat
  }

  # Q3 for observed data -------------------------------------------------------
  q3 <- q3_func(object$y, p)

  # Q3 for replicated datasets -------------------------------------------------
  rep_q3 <- vector("list", length = total_draws)
  for (i in 1:total_draws) {
    rep_q3[[i]] <- q3_func(rep_y[[i]], rep_p[[i]])
  }

  # Compute Q3 PPP -------------------------------------------------------------
  ppp_q3 <- matrix(NA, ncol(object$y), ncol(object$y))
  for (i in 1:ncol(object$y)) {
    for (ii in 1:ncol(object$y)) {
      ppp_q3[i, ii] <- mean(x = vapply(
        X = rep_q3,
        FUN = function(m) m[i, ii], numeric(1)
      ) > q3[i, ii], na.rm = TRUE)
    }
  }

  # Return results of function -------------------------------------------------
  list(
    q3 = q3,
    ppp_q3 = ppp_q3
  )
}
