#-------------------------------------------------------------------------------
#' MCMH Parameter Estimates for Multiple Chains
#'
#' This function calculates MCMH parameter estimates for multiple chains. See
#' documentation for mcmhSC.R for more information.
#'
#' @examples
#'mcmhMC(chains = 3, y = sdirt$y, objFun = dich_response_model, estOmega = T,
#'     estNu = T, estZeta = T, lambda = sdirt$lambda, gamma = sdirt$gamma,
#'     omega0 = array(data = 0, dim = dim(sdirt$omega)),
#'     nu0 = array(data = 0, dim = c(ncol(sdirt$nu), 1)),
#'     zeta0 = array(data = 0, dim = dim(sdirt$zeta)),
#'     omega_mu = sdirt$omega_mu, omega_sigma2 = sdirt$omega_sigma2,
#'     nu_mu = matrix(sdirt$nu_mu), nu_sigma2 = matrix(sdirt$nu_sigma2),
#'     zeta_mu = sdirt$zeta_mu, zeta_sigma2 = sdirt$zeta_sigma2,
#'     burn=0, thin=10, min_tune=50, tune_int=50, max_tune = 1000, niter = 2000,
#'     weight=1/1)
#'
#' @export mcmhMC
#-------------------------------------------------------------------------------

mcmhMC <- function(
  chains=NULL, y = y, objFun = NULL, estOmega = T, estNu = T, estZeta = T,
  lambda = NULL, gamma = NULL, omega0 = NULL, nu0 = NULL, zeta0 = NULL,
  omega_mu = NULL, omega_sigma2 = NULL, nu_mu = NULL, nu_sigma2 = NULL,
  zeta_mu = NULL, zeta_sigma2 = NULL, burn = NULL, thin = NULL, min_tune = NULL,
  tune_int = NULL, max_tune = NULL, niter = NULL, weight = NULL
)
{
  draws = parallel::mclapply(
    mc.cores = parallel::detectCores(), X = 1:chains, FUN = mcmhSC, y = y,
    objFun = objFun, estOmega = estOmega, estNu = estNu, estZeta = estZeta,
    lambda = lambda, gamma = gamma, omega0 = omega0, nu0 = nu0, zeta0 = zeta0,
    omega_mu = omega_mu, omega_sigma2 = omega_sigma2, nu_mu = nu_mu,
    nu_sigma2 = nu_sigma2, zeta_mu = zeta_mu, zeta_sigma2 = zeta_sigma2,
    burn = burn, thin = thin, min_tune = min_tune, tune_int = tune_int,
    max_tune = max_tune, niter = niter, weight = weight
  )
  names(draws) <- paste("chain",1:chains,sep = "")
  if(estOmega) {
    omegadraws <- abind::abind(
      lapply(X = draws, FUN = function(x) {x[["omega_draws"]]}),
      along = 1
    )
    omegaeap <- t(apply(X = omegadraws, MARGIN = c(2,3), FUN = mean))
    omegapsd <- t(apply(X = omegadraws, MARGIN = c(2,3), FUN = sd))
    if(chains > 1) {
      omegaPSRF <- matrix(
        data = unlist(
          lapply(X = 1:nrow(omega0), FUN = function(i) {
            coda::gelman.diag(
              x = coda::mcmc.list(
                lapply(
                  X = 1:chains,
                  FUN = function(ch) {
                    coda::mcmc(data = lapply(
                      X = draws,
                      FUN = function(x) {x[["omega_draws"]]}
                    )[[ch]][,,i])
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
      omegaPSRF <- NULL
    }
  } else {
    omegaeap <- NULL
    omegapsd <- NULL
    omegaPSRF <- NULL
  }
  if(estNu) {
    nudraws <- abind::abind(
      lapply(X = draws, FUN = function(x) {x[["nu_draws"]]}),
      along = 1
    )
    nueap <- t(apply(X = nudraws, MARGIN = c(2,3), FUN = mean))
    nupsd <- t(apply(X = nudraws, MARGIN = c(2,3), FUN = sd))
    if(chains > 1) {
      nuPSRF <- matrix(
        data = unlist(
          lapply(X = 1:nrow(nu0), FUN = function(i) {
            coda::gelman.diag(
              x = coda::mcmc.list(
                lapply(
                  X = 1:chains,
                  FUN = function(ch) {
                    coda::mcmc(data = lapply(
                      X = draws,
                      FUN = function(x) {x[["nu_draws"]]}
                    )[[ch]][,,i])
                  }
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
      nuPSRF <- NULL
    }
  } else {
    nueap <- NULL
    nupsd <- NULL
    nuPSRF <- NULL
  }
  if(estZeta) {
    zetadraws <- abind::abind(
      lapply(X = draws, FUN = function(x) {x[["zeta_draws"]]}),
      along = 1
    )
    zetaeap <- t(apply(X = zetadraws, MARGIN = c(2,3), FUN = mean))
    zetapsd <- t(apply(X = zetadraws, MARGIN = c(2,3), FUN = sd))
    if(chains>1) {
      zetaPSRF <- matrix(
        data = unlist(
          lapply(X = 1:nrow(zeta0), FUN = function(i) {
            coda::gelman.diag(
              x = coda::mcmc.list(
                lapply(
                  X = 1:chains,
                  FUN = function(ch) {
                    coda::mcmc(data = lapply(
                      X = draws,
                      FUN = function(x) {x[["zeta_draws"]]}
                    )[[ch]][,,i])
                  }
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
      zetaPSRF <- NULL
    }
  } else {
    zetaeap <- NULL
    zetapsd <- NULL
    zetaPSRF <- NULL
  }
  return(list(
    "mcmhdraws" = draws, "omegaEAP" = omegaeap, "omegaPSD" = omegapsd,
    "omegaPSRF" = omegaPSRF, "nuEAP" = nueap, "nuPSD" = nupsd,
    "nuPSRF" = nuPSRF, "zetaEAP" = zetaeap, "zetaPSD" = zetapsd,
    "zetaPSRF" = zetaPSRF
  ))
}
