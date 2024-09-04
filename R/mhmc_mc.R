#-------------------------------------------------------------------------------
#' MHMC Parameter Estimates for Multiple Chains
#'
#' This function calculates MHMC parameter estimates for multiple chains. See
#' documentation for mhmc_sc.R for more information.
#'
#' @param chains Number of chains in the MHMC sampler (scalar).
#' @param y Matrix of item responses (K by IJ).
#' @param obj_fun A function that calculates predictions and log-likelihood
#' values for the selected model (character).
#' @param link Choose between "logit" or "probit" link functions.
#' @param est_omega Determines whether omega is estimated (logical).
#' @param est_lambda Determines whether nu is estimated (logical).
#' @param est_zeta Determines whether zeta is estimated (logical).
#' @param est_nu Determines whether nu is estimated (logical).
#' @param omega0 Starting or known values for omega (K by MN).
#' @param gamma0 Starting or known values for gamma (JM by MN).
#' @param lambda0 Starting or known values for lambda (IJ by JM).
#' @param zeta0 Starting or known values for  zeta (K by JM).
#' @param nu0 Starting or known values for nu (IJ by 1).
#' @param kappa0 Starting or known values for kappa (1 by IJ).
#' @param omega_mu Mean prior for omega (1 by MN).
#' @param omega_sigma2 Covariance prior for omega (MN by MN).
#' @param lambda_mu Mean prior for lambda (1 by JM)
#' @param lambda_sigma2 Covariance prior for lambda (JM by JM)
#' @param zeta_mu Mean prior for zeta (1 by JM).
#' @param zeta_sigma@ Covariance prior for zeta (JM by JM).
#' @param nu_mu Mean prior for nu (1 by 1)
#' @param nu_sigma2 Covariance prior for nu (1 by 1)
#' @param burn Number of iterations at the beginning of an MCMC run to discard
#' (scalar).
#' @param thin Determines every nth observation retained (scalar).
#' @param min_tune Determines when tunning begins (scalar).
#' @param tune_int MHMC tuning interval (scalar).
#' @param max_tune Determines when tunning ends (scalar).
#' @param niter Number of iterations of the MHMC sampler.
#' @param psrf Estimate potential scale reduction factor (logical).
#'
#' @return List with elements omega_draws (draws from every saved iteration of
#' the MHMC sampler), omegaEAP (expected a posteriori estimates for omega),
#' omegaPSD (posterior standard deviation estimates for omega), omega_psrf
#' (potential scale reduction factor for omega), nuEAP (expected a posteriori
#' estimates for nu), nuPSD (posterior standard deviation estimates for nu),
#' nu_psrf (potential scale reduction factor for nu), zetaEAP (expected a
#' posteriori estimates for zeta), zetaPSD (posterior standard deviation
#' estimates for zeta), zeta_psrf (potential scale reduction factor for zeta).
#'
#' @examples
#' mhmc_mc(chains = 3, y = ex1$y, obj_fun = dich_response_model,
#'         est_omega = TRUE, est_lambda = FALSE, est_zeta = TRUE, est_nu = TRUE,
#'         omega0 = array(data = 0, dim = dim(ex1$omega)), gamma0 = ex1$gamma,
#'         lambda0 = ex1$lambda, zeta0 = array(data = 0, dim = dim(ex1$zeta)),
#'         nu0 = array(data = 0, dim = c(ncol(ex1$y), 1)),
#'         omega_mu = ex1$omega_mu, omega_sigma2 = ex1$omega_sigma2,
#'         zeta_mu = ex1$zeta_mu, zeta_sigma2 = ex1$zeta_sigma2,
#'         nu_mu = ex1$nu_mu, nu_sigma2 = ex1$nu_sigma2, burn = 0, thin = 10,
#'         min_tune = 50, tune_int = 50, max_tune = 1000, niter = 2000,
#'         psrf = TRUE)
#'
#' @keywords internal
#-------------------------------------------------------------------------------

mhmc_mc <- function(
    chains = NULL,
    y = y,
    obj_fun = NULL,
    link = NULL,
    est_omega = TRUE,
    est_lambda = TRUE,
    est_zeta = TRUE,
    est_nu = TRUE,
    omega0 = NULL,
    gamma0 = NULL,
    lambda0 = NULL,
    zeta0 = NULL,
    nu0 = NULL,
    kappa0 = NULL,
    omega_mu = NULL,
    omega_sigma2 = NULL,
    lambda_mu = NULL,
    lambda_sigma2 = NULL,
    zeta_mu = NULL,
    zeta_sigma2 = NULL,
    nu_mu = NULL,
    nu_sigma2 = NULL,
    burn = NULL,
    thin = NULL,
    min_tune = NULL,
    tune_int = NULL,
    max_tune = NULL,
    niter = NULL,
    psrf = FALSE
) {
  if (!requireNamespace("abind", quietly = TRUE)) {
    stop("Package \"abind\" needed for the mhmc_mc function to work. Please
         install.",
         call. = FALSE)
    stop("Package \"parallel\" needed for the mhmc_mc function to work. Please
         install.",
         call. = FALSE)
    stop("Package \"coda\" needed for the mhmc_mc function to work. Please
         install.",
         call. = FALSE)
  }
  draws <- parallel::mclapply(
    mc.cores = parallel::detectCores(), X = 1:chains, FUN = mhmc_sc, y = y,
    obj_fun = obj_fun, link = link, est_omega = est_omega,
    est_lambda = est_lambda, est_zeta = est_zeta, est_nu = est_nu,
    omega0 = omega0, gamma0 = gamma0, lambda0 = lambda0, zeta0 = zeta0,
    nu0 = nu0, kappa0 = kappa0, omega_mu = omega_mu,
    omega_sigma2 = omega_sigma2, lambda_mu = lambda_mu,
    lambda_sigma2 = lambda_sigma2, zeta_mu = zeta_mu, zeta_sigma2 = zeta_sigma2,
    nu_mu = nu_mu, nu_sigma2 = nu_sigma2, burn = burn,
    thin = thin, min_tune = min_tune, tune_int = tune_int, max_tune = max_tune,
    niter = niter
  )
  names(draws) <- paste("chain", 1:chains, sep = "")
  if (est_omega) {
    omegadraws <- abind::abind(
      lapply(X = draws, FUN = function(x) {
        x[["omega_draws"]]
      }
      ),
      along = 1
    )
    omegaeap <- t(apply(X = omegadraws, MARGIN = c(2, 3), FUN = mean))
    omegapsd <- t(apply(X = omegadraws, MARGIN = c(2, 3), FUN = sd))
    if (chains > 1 && psrf) {
      omega_psrf <- matrix(
        data = unlist(
          lapply(X = seq_len(length.out = nrow(omega0)), FUN = function(i) {
            coda::gelman.diag(
              x = coda::mcmc.list(
                lapply(
                  X = 1:chains,
                  FUN = function(ch) {
                    coda::mcmc(data = lapply(
                      X = draws,
                      FUN = function(x) {
                        x[["omega_draws"]]
                      }
                    )[[ch]][, , i])
                  }
                )
              ),
              multivariate = FALSE
            )[["psrf"]][, 1]})
        ),
        nrow = nrow(omega0),
        ncol = ncol(omega0),
        byrow = TRUE
      )
    } else {
      omega_psrf <- NULL
    }
  } else {
    omegaeap <- NULL
    omegapsd <- NULL
    omega_psrf <- NULL
  }
  if (est_lambda) {
    lambdadraws <- abind::abind(
      lapply(X = draws, FUN = function(x) {
        x[["lambda_draws"]]
      }
      ),
      along = 1
    )
    lambdaeap <- t(apply(X = lambdadraws, MARGIN = c(2, 3), FUN = mean))
    lambdapsd <- t(apply(X = lambdadraws, MARGIN = c(2, 3), FUN = sd))
    if (chains > 1 && psrf) {
      lambda_psrf <- matrix(
        data = unlist(
          lapply(X = seq_len(length.out = nrow(lambda0)), FUN = function(i) {
            coda::gelman.diag(
              x = coda::mcmc.list(
                lapply(
                  X = 1:chains,
                  FUN = function(ch) {
                    coda::mcmc(data = lapply(
                      X = draws,
                      FUN = function(x) {
                        x[["lambda_draws"]]
                      }
                    )[[ch]][, , i])
                  }
                )
              ),
              multivariate = FALSE
            )[["psrf"]][, 1]})
        ),
        nrow = nrow(x = lambda0),
        ncol = ncol(x = lambda0),
        byrow = TRUE
      )
    } else {
      lambda_psrf <- NULL
    }
  } else {
    lambdaeap <- NULL
    lambdapsd <- NULL
    lambda_psrf <- NULL
  }
  if (est_zeta) {
    zetadraws <- abind::abind(
      lapply(X = draws, FUN = function(x) {
        x[["zeta_draws"]]
      }),
      along = 1
    )
    zetaeap <- t(apply(X = zetadraws, MARGIN = c(2, 3), FUN = mean))
    zetapsd <- t(apply(X = zetadraws, MARGIN = c(2, 3), FUN = sd))
    if (chains > 1 && psrf) {
      zeta_psrf <- matrix(
        data = unlist(
          lapply(X = seq_len(nrow(zeta0)), FUN = function(i) {
            coda::gelman.diag(
              x = coda::mcmc.list(
                lapply(
                  X = 1:chains,
                  FUN = function(ch) {
                    coda::mcmc(
                      data = lapply(X = draws, FUN = function(x) {
                        x[["zeta_draws"]]
                      })[[ch]][, , i]
                    )
                  }
                )
              ),
              multivariate = FALSE
            )[["psrf"]][, 1]})
        ),
        nrow = nrow(x = zeta0),
        ncol = ncol(x = zeta0),
        byrow = TRUE
      )
    } else {
      zeta_psrf <- NULL
    }
  } else {
    zetaeap <- NULL
    zetapsd <- NULL
    zeta_psrf <- NULL
  }
  if (est_nu) {
    nudraws <- abind::abind(
      lapply(X = draws, FUN = function(x) {
        x[["nu_draws"]]
      }),
      along = 1
    )
    nueap <- t(apply(X = nudraws, MARGIN = c(2, 3), FUN = mean))
    nupsd <- t(apply(X = nudraws, MARGIN = c(2, 3), FUN = sd))
    if (chains > 1 && psrf) {
      nu_psrf <- matrix(
        data = unlist(
          lapply(X = seq_len(nrow(nu0)), FUN = function(i) {
            coda::gelman.diag(
              x = coda::mcmc.list(
                lapply(
                  X = 1:chains,
                  FUN = function(ch) {
                    coda::mcmc(data = lapply(X = draws, FUN = function(x) {
                      x[["nu_draws"]]
                    })[[ch]][, , i]
                    )
                  }
                )
              ),
              multivariate = FALSE
            )[["psrf"]][, 1]})
        ),
        nrow = nrow(nu0),
        ncol = 1,
        byrow = TRUE
      )
    } else {
      nu_psrf <- NULL
    }
  } else {
    nueap <- NULL
    nupsd <- NULL
    nu_psrf <- NULL
  }
  return(list(
    "mhmcdraws" = draws, "omegaEAP" = omegaeap, "omegaPSD" = omegapsd,
    "omega_psrf" = omega_psrf,  "lambdaEAP" = lambdaeap,
    "lambdaPSD" = lambdapsd, "lambda_psrf" = lambda_psrf,  "zetaEAP" = zetaeap,
    "zetaPSD" = zetapsd, "zeta_psrf" = zeta_psrf, "nuEAP" = nueap,
    "nuPSD" = nupsd, "nu_psrf" = nu_psrf
  ))
}
