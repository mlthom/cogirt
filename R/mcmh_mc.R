#-------------------------------------------------------------------------------
#' MCMH Parameter Estimates for Multiple Chains
#'
#' This function calculates MCMH parameter estimates for multiple chains. See
#' documentation for mcmh_sc.R for more information.
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
#'
#' @return List with elements omega_draws (draws from every saved iteration of
#' the MCMH sampler), omegaEAP (expected a posteriori estimates for omega),
#' omegaPSD (posterior standard deviation estimates for omega), omega_psrf
#' (potential scale reduction factor for omega), nuEAP (expected a posteriori
#' estimates for nu), nuPSD (posterior standard deviation estimates for nu),
#' nu_psrf (potential scale reduction factor for nu), zetaEAP (expected a
#' posteriori estimates for zeta), zetaPSD (posterior standard deviation
#' estimates for zeta), zeta_psrf (potential scale reduction factor for zeta).
#'
#' @examples
#'mcmh_mc(chains = 3, y = sdirt$y, obj_fun = dich_response_model, est_omega = T,
#'     est_nu = T, est_zeta = T, lambda = sdirt$lambda, kappa = sdirt$kappa,
#'     gamma = sdirt$gamma, omega0 = array(data = 0, dim = dim(sdirt$omega)),
#'     nu0 = array(data = 0, dim = c(ncol(sdirt$nu), 1)),
#'     zeta0 = array(data = 0, dim = dim(sdirt$zeta)),
#'     omega_mu = sdirt$omega_mu, omega_sigma2 = sdirt$omega_sigma2,
#'     nu_mu = matrix(sdirt$nu_mu), nu_sigma2 = matrix(sdirt$nu_sigma2),
#'     zeta_mu = sdirt$zeta_mu, zeta_sigma2 = sdirt$zeta_sigma2,
#'     burn = 0, thin = 10, min_tune = 50, tune_int = 50, max_tune = 1000,
#'     niter = 2000)
#'
#' @export mcmh_mc
#-------------------------------------------------------------------------------

mcmh_mc <- function(
  chains=NULL, y = y, obj_fun = NULL, est_omega = T, est_nu = T, est_zeta = T,
  lambda = NULL, kappa = NULL, gamma = NULL, omega0 = NULL, nu0 = NULL,
  zeta0 = NULL, omega_mu = NULL, omega_sigma2 = NULL, nu_mu = NULL,
  nu_sigma2 = NULL, zeta_mu = NULL, zeta_sigma2 = NULL, burn = NULL,
  thin = NULL, min_tune = NULL, tune_int = NULL, max_tune = NULL, niter = NULL
) {
  if (!requireNamespace("abind", quietly = TRUE)) {
    stop("Package \"abind\" needed for the mcmh_mc function to work. Please
         install.",
         call. = FALSE)
    stop("Package \"parallel\" needed for the mcmh_mc function to work. Please
         install.",
         call. = FALSE)
    stop("Package \"coda\" needed for the mcmh_mc function to work. Please
         install.",
         call. = FALSE)
  }
  draws <- parallel::mclapply(
    mc.cores = parallel::detectCores(), X = 1:chains, FUN = mcmh_sc, y = y,
    obj_fun = obj_fun, est_omega = est_omega, est_nu = est_nu,
    est_zeta = est_zeta, lambda = lambda, kappa = kappa, gamma = gamma,
    omega0 = omega0, nu0 = nu0, zeta0 = zeta0, omega_mu = omega_mu,
    omega_sigma2 = omega_sigma2, nu_mu = nu_mu, nu_sigma2 = nu_sigma2,
    zeta_mu = zeta_mu, zeta_sigma2 = zeta_sigma2, burn = burn, thin = thin,
    min_tune = min_tune, tune_int = tune_int, max_tune = max_tune, niter = niter
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
    if (chains > 1) {
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
              multivariate = F
            )[["psrf"]][, 1]})
        ),
        nrow = nrow(omega0),
        ncol = ncol(omega0),
        byrow = T
      )
    } else {
      omega_psrf <- NULL
    }
  } else {
    omegaeap <- NULL
    omegapsd <- NULL
    omega_psrf <- NULL
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
    if (chains > 1) {
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
                    )}
                )
              ),
              multivariate = F
            )[["psrf"]][, 1]})
        ),
        nrow = nrow(nu0),
        ncol = 1,
        byrow = T
      )
    } else {
      nu_psrf <- NULL
    }
  } else {
    nueap <- NULL
    nupsd <- NULL
    nu_psrf <- NULL
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
    if (chains > 1) {
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
                    )}
                )
              ),
              multivariate = F
            )[["psrf"]][, 1]})
        ),
        nrow = nrow(zeta0),
        ncol = ncol(zeta0),
        byrow = T
      )
    } else {
      zeta_psrf <- NULL
    }
  } else {
    zetaeap <- NULL
    zetapsd <- NULL
    zeta_psrf <- NULL
  }
  return(list(
    "mcmhdraws" = draws, "omegaEAP" = omegaeap, "omegaPSD" = omegapsd,
    "omega_psrf" = omega_psrf, "nuEAP" = nueap, "nuPSD" = nupsd,
    "nu_psrf" = nu_psrf, "zetaEAP" = zetaeap, "zetaPSD" = zetapsd,
    "zeta_psrf" = zeta_psrf
  ))
}
