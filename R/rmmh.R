#-------------------------------------------------------------------------------
#' RMMH Parameter Estimates for Multiple Chains
#'
#' This function calculates RMMH parameter estimates for multiple chains.
#'
#' @param chains Number of chains in the MCMH sampler (scalar).
#' @param y Matrix of item responses (K by IJ).
#' @param obj_fun A function that calculates predictions and log-likelihood
#' values for the selected model (character).
#' @param est_omega Determines whether omega is estimated (logical).
#' @param est_nu Determines whether nu is estimated (logical).
#' @param est_zeta Determines whether zeta is estimated (logical).
#' @param lambda Matrix of item structure parameters (IJ by JM).
#' @param kappa Matrix of item guessing parameters (K by IJ).
#' @param gamma Matrix of experimental structure parameters (JM by MN).
#' @param omega0 Starting values for omega.
#' @param nu0 Starting values for nu.
#' @param zeta0 Starting values for zeta.
#' @param omega_mu Vector of means prior for omega (1 by MN).
#' @param omega_sigma Covariance matrix prior for omega (MN by MN).
#' @param zeta_mu Vector of means prior for zeta (1 by JM).
#' @param zeta_sigma@ Covariance matrix prior for zeta (JM by JM).
#' @param nu_mu Prior mean for nu (scalar).
#' @param nu_sigma@ Prior variance for nu (scalar).
#' @param burn Number of iterations at the beginning of an MCMC run to discard
#' (scalar).
#' @param thin Determines every nth observation retained (scalar).
#' @param min_tune Determines when tunning begins (scalar).
#' @param tune_int MCMH tuning interval (scalar).
#' @param max_tune Determines when tunning ends (scalar).
#' @param niter Number of iterations of the MCMH sampler.
#' @param verbose_rmmh Print progress of MCMH sampler.
#' @param max_iter_rmmh Maximum number of iterations for RMMH.
#'
#' @references
#' Cai, L. (2010). High-dimensional exploratory item factor analysis by a
#' Metropolis-Hastings Robbins-Monro algorithm. \emph{Psychometrika, 75(1)},
#' 33-57.
#'
#' Cai, L. (2010). Metropolis-Hastings Robbins-Monro algorithm for confirmatory
#' item factor analysis. \emph{Journal of Educational and Behavioral Statistics,
#' 35(3)}, 307-335.
#'
#' @examples
#' # Multiple subjects example. Extended MCMH sampling.
#' rmmh(chains = 3, y = sdirt$y, obj_fun = dich_response_model, est_omega = T,
#'     est_nu = T, est_zeta = T, lambda = sdirt$lambda, kappa = sdirt$kappa,
#'     gamma = sdirt$gamma, omega0 = array(data = 0, dim = dim(sdirt$omega)),
#'     nu0 = array(data = 0, dim = c(ncol(sdirt$nu), 1)),
#'     zeta0 = array(data = 0, dim = dim(sdirt$zeta)),
#'     omega_mu = sdirt$omega_mu, omega_sigma2 = sdirt$omega_sigma2,
#'     nu_mu = matrix(sdirt$nu_mu), nu_sigma2 = matrix(sdirt$nu_sigma2),
#'     zeta_mu = sdirt$zeta_mu, zeta_sigma2 = sdirt$zeta_sigma2,
#'     burn=50, thin=10, min_tune=10, tune_int=10, max_tune = 100,
#'     niter = 100, verbose_rmmh = T, max_iter_rmmh = 200)
#'
#' # Multiple subjects example. Limited MCMH sampling.
#' rmmh(chains = 3, y = sdirt$y, obj_fun = dich_response_model, est_omega = T,
#'     est_nu = T, est_zeta = T, lambda = sdirt$lambda, kappa = sdirt$kappa,
#'     gamma = sdirt$gamma, omega0 = array(data = 0, dim = dim(sdirt$omega)),
#'     nu0 = array(data = 0, dim = c(ncol(sdirt$nu), 1)),
#'     zeta0 = array(data = 0, dim = dim(sdirt$zeta)),
#'     omega_mu = sdirt$omega_mu, omega_sigma2 = sdirt$omega_sigma2,
#'     nu_mu = matrix(sdirt$nu_mu), nu_sigma2 = matrix(sdirt$nu_sigma2),
#'     zeta_mu = sdirt$zeta_mu, zeta_sigma2 = sdirt$zeta_sigma2,
#'     burn=0, thin=1, min_tune=0, tune_int=0, max_tune = 0, niter = 1,
#'     verbose_rmmh = T, max_iter_rmmh = 200)
#'
#' # Single subject example. Extended MCMH sampling.
#' rmmh(chains = 3, y = sdirtSS$y, obj_fun = dich_response_model, est_omega = T,
#'      est_nu = T, est_zeta = T, lambda = sdirtSS$lambda, kappa = sdirt$kappa,
#'      gamma = sdirtSS$gamma, omega0 = array(data = 0,
#'      dim = dim(sdirtSS$omega)), nu0 = array(data = 0,
#'      dim = c(ncol(sdirtSS$nu), 1)),
#'      zeta0 = array(data = 0, dim = dim(sdirtSS$zeta)),
#'      omega_mu = sdirtSS$omega_mu, omega_sigma2 = sdirtSS$omega_sigma2,
#'      nu_mu = matrix(sdirtSS$nu_mu), nu_sigma2 = matrix(sdirtSS$nu_sigma2),
#'      zeta_mu = sdirtSS$zeta_mu, zeta_sigma2 = sdirtSS$zeta_sigma2,
#'      burn=50, thin=10, min_tune=10, tune_int=10, max_tune = 100,
#'      niter = 100, verbose_rmmh = T, max_iter_rmmh = 200)
#'
#' # Single subject example. Limited MCMH sampling.
#' rmmh(chains = 3, y = sdirtSS$y, obj_fun = dich_response_model,
#'      est_omega = T, est_nu = T, est_zeta = T, lambda = sdirtSS$lambda,
#'      kappa = sdirt$kappa, gamma = sdirtSS$gamma,
#'      omega0 = array(data = 0, dim = dim(sdirtSS$omega)),
#'      nu0 = array(data = 0, dim = c(ncol(sdirtSS$nu), 1)),
#'      zeta0 = array(data = 0, dim = dim(sdirtSS$zeta)),
#'      omega_mu = sdirtSS$omega_mu, omega_sigma2 = sdirtSS$omega_sigma2,
#'      nu_mu = matrix(sdirtSS$nu_mu), nu_sigma2 = matrix(sdirtSS$nu_sigma2),
#'      zeta_mu = sdirtSS$zeta_mu, zeta_sigma2 = sdirtSS$zeta_sigma2,
#'      burn=0, thin=1, min_tune=0, tune_int=0, max_tune = 0, niter = 1,
#'      verbose_rmmh = T, max_iter_rmmh = 200)
#'
#' @export rmmh
#-------------------------------------------------------------------------------

rmmh <- function(
  chains = 1, y = y, obj_fun = NULL, est_omega = T, est_nu = T, est_zeta = T,
  lambda = NULL, kappa = NULL, gamma = NULL, omega0 = NULL, nu0 = NULL,
  zeta0 = NULL, omega_mu = NULL, omega_sigma2 = NULL, nu_mu = NULL,
  nu_sigma2 = NULL, zeta_mu = NULL, zeta_sigma2 = NULL, burn = NULL,
  thin = NULL, min_tune = NULL, tune_int = NULL, max_tune = NULL,
  niter = NULL, verbose_rmmh = T, max_iter_rmmh = 200
) {
  if (!requireNamespace("abind", quietly = TRUE)) {
    stop("Package \"abind\" needed for the rmmh function to work. Please
         install.",
         call. = FALSE)
  }
  # STEP 0: MCMC burn in -------------------------------------------------------
  if (verbose_rmmh) {
    cat(
      "MCMC Burn-In Start Time",
      format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
      "\n",
      sep = " "
    )
  }
  mcmhburn <- mcmh_mc(
    chains = chains, y = y, obj_fun = obj_fun, est_omega = est_omega,
    est_nu = est_nu, est_zeta = est_zeta, lambda = lambda, kappa = kappa,
    gamma = gamma, omega0 = omega0, nu0 = nu0, zeta0 = zeta0,
    omega_mu = omega_mu, omega_sigma2 = omega_sigma2, nu_mu = nu_mu,
    nu_sigma2 = nu_sigma2, zeta_mu = zeta_mu, zeta_sigma2 = zeta_sigma2,
    burn = burn, thin = thin, min_tune = min_tune, tune_int = tune_int,
    max_tune = max_tune, niter = niter
  )

  # Update initial estimates and variance of candidates
  if (est_omega) {
    omega0 <- mcmhburn$omegaEAP
  }
  if (est_nu) {
    nu0 <- mcmhburn$nuEAP
  }
  if (est_zeta) {
    zeta0 <- mcmhburn$zetaEAP
  }
  # Set up RMMH loop parameters
  tol <- .01
  log_lik <- 0
  iter <- 0
  info0 <- diag(x = 10, nrow = ncol(omega0))
  go <- T
  ongoing_omega <- omega0
  #start iteration
  if (verbose_rmmh) {
    cat(
      "... burn-in completed at",
      format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
      "\n",
      sep = " "
    )
    cat(
      "RM-MH Start Time",
      format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
      "\n",
      sep = " "
    )
  }

  while (go && (iter < max_iter_rmmh)) {
    iter <- iter + 1

    # STEP 1: Stochastic imputation --------------------------------------------
    mc_draws_at_iteration_k <- mcmh_mc(
      chains = 3, y = y, obj_fun = obj_fun, est_omega = est_omega,
      est_nu = est_nu, est_zeta = est_zeta, lambda = lambda, kappa = kappa,
      gamma = gamma, omega0 = omega0, nu0 = nu0, zeta0 = zeta0,
      omega_mu = omega_mu, omega_sigma2 = omega_sigma2, nu_mu = matrix(nu_mu),
      nu_sigma2 = matrix(nu_sigma2), zeta_mu = zeta_mu,
      zeta_sigma2 = zeta_sigma2, burn = burn, thin = thin, min_tune = min_tune,
      tune_int = tune_int, max_tune = max_tune, niter = niter
    )

    # STEP 2: Stochastic approximation -----------------------------------------
    grad <- NULL
    for (i in 1:(length(seq(from = burn + 1,
                            to = niter,
                            by = thin)) * chains)) {
      grad <-
        abind::abind(
          grad,
          matrix(
            data = unlist(
              dich_response_deriv(
                y = y,
                nu = array(
                  data = if (est_nu) {
                    array(
                      data = abind::abind(
                        lapply(mc_draws_at_iteration_k$mcmhdraws, function(x) {
                          x[["nu_draws"]]
                        }),
                        along = 1
                      )[i, , ],
                      dim = dim(nu0)
                    )
                  } else {
                    nu0
                  },
                  dim = dim(y)
                ),
                lambda = lambda,
                kappa = kappa,
                gamma = gamma,
                omega = omega0,
                zeta = if (est_zeta) {
                  matrix(
                    data = abind::abind(
                      lapply(mc_draws_at_iteration_k$mcmhdraws, function(x) {
                        x[["zeta_draws"]]
                      }),
                      along = 1
                    )[i, , ],
                    nrow = nrow(zeta0),
                    ncol = ncol(zeta0),
                    byrow = T
                  )
                } else{
                  zeta0
                },
                omega_mu = omega_mu,
                omega_sigma2 = omega_sigma2,
                zeta_mu = zeta_mu,
                zeta_sigma2 = zeta_sigma2
              )$fpd),
            nrow = nrow(omega0),
            ncol = ncol(omega0),
            byrow = T
          ),
          along = 3
        )
    }
    grad <- apply(X = grad, MARGIN = c(1, 2), FUN = mean)
    info <- NULL
    for (i in 1:(length(seq(from = burn + 1,
                            to = niter,
                            by = thin)) * chains)) {
      info <- abind::abind(
        info,
        array(
          data = unlist(
            dich_response_deriv(
              y = y,
              nu = array(
                data = if (est_nu) {
                  array(
                    data = abind::abind(
                      lapply(
                        X = mc_draws_at_iteration_k$mcmhdraws,
                        FUN = function(x) {
                          x[["nu_draws"]]
                        }),
                      along = 1)[i, , ],
                    dim = dim(nu0)
                  )
                } else {
                  nu0
                },
                dim = dim(y)
              ),
              lambda = lambda,
              kappa = kappa,
              gamma = gamma,
              omega = omega0,
              zeta = if (est_zeta) {
                matrix(
                  data = abind::abind(
                    lapply(
                      X = mc_draws_at_iteration_k$mcmhdraws,
                      FUN = function(x) {
                        x[["zeta_draws"]]
                      }),
                    along = 1)[i, , ],
                  nrow = nrow(zeta0),
                  ncol = ncol(zeta0),
                  byrow = T)
              } else{
                zeta0
              },
              omega_mu = omega_mu,
              omega_sigma2 = omega_sigma2,
              zeta_mu = zeta_mu,
              zeta_sigma2 = zeta_sigma2
            )$post_info), c(ncol(omega0), ncol(omega0), nrow(y))),
        along = 3
      )
    }
    info <- apply(X = info, MARGIN = c(1, 2), FUN = mean)
    gain <- 1 / iter
    info1 <- info0 + gain * (info - info0)

    # STEP 3: Robbins-Monro update ---------------------------------------------
    # Note using ginv instead of solve because appears to be more stable. Later
    # added error handling for cases when matrix is singular.
    inv <- try(expr = MASS::ginv(X = info1, tol = .Machine$double.xmin))
    if (all(class(x = inv) == "try-error")) {
      omega1 <- mcmhburn$omegaEAP
    } else {
      omega1 <- omega0 + gain * grad %*% inv
    }
    p <- obj_fun(y = y, nu = array(data = nu0, dim = dim(y)), lambda = lambda,
                 kappa = kappa, gamma = gamma, omega = omega1, zeta = zeta0)$p
    log_lik <- sum(x = log((p^y) * (1 - p) ^ (1 - y)), na.rm = T)

    # Test for completion
    if (any(abs(omega1 - omega0) > tol) & iter < max_iter_rmmh) {
      go <- T
      #update values
      info0 <- info1
      omega0 <- omega1
      ongoing_omega <- c(ongoing_omega, omega1)
      if (est_nu) {
        nu0 <- mc_draws_at_iteration_k$nuEAP
      }
      if (est_zeta) {
        zeta0 <- mc_draws_at_iteration_k$zetaEAP
      }
      if (verbose_rmmh) {
        cat(
          "\r",
          "... at iteration ",
          iter,
          " MAP omega is ",
          format(x = round(x = omega1, digits = 3), nsmall = 3),
          " logLik is ",
          format(x = round(x = log_lik, digits = 4), nsmall = 4),
          sep = ""
        )
      }
    } else if (all(abs(omega1 - omega0) < tol) & iter < max_iter_rmmh) {
      go <- F
      if (verbose_rmmh) {
        cat(
          "\n... algorithm converged at",
          format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
          "\n",
          sep = " "
        )
      }
    } else if (all(abs(omega1 - omega0) < tol) & iter == max_iter_rmmh) {
      if (verbose_rmmh) {
        cat(
          "\n... algorithm failed to converge at",
          format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
          "\n",
          sep = " "
        )
      }
    }
  }
  return(list(
    "omega1" = omega1,
    "info1" = info1,
    "log_lik" = log_lik,
    "ongoing_omega" = ongoing_omega
  ))
}
