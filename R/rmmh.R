#-------------------------------------------------------------------------------
#' RMMH Parameter Estimates for Multiple Chains
#'
#' This function calculates RMMH parameter estimates for multiple chains.
#'
#' @references

#'
#' @examples
#'rmmh(chains = 3, y = sdirt$y, objFun = dich_response_model, estOmega = T,
#'     estNu = T, estZeta = T, lambda = sdirt$lambda, gamma = sdirt$gamma,
#'     omega0 = array(data = 0, dim = dim(sdirt$omega)),
#'     nu0 = array(data = 0, dim = c(ncol(sdirt$nu), 1)),
#'     zeta0 = array(data = 0, dim = dim(sdirt$zeta)),
#'     omega_mu = sdirt$omega_mu, omega_sigma2 = sdirt$omega_sigma2,
#'     nu_mu = matrix(sdirt$nu_mu), nu_sigma2 = matrix(sdirt$nu_sigma2),
#'     zeta_mu = sdirt$zeta_mu, zeta_sigma2 = sdirt$zeta_sigma2,
#'     burn=0, thin=10, min_tune=50, tune_int=50, max_tune = 1000, niter = 2000,
#'     weight=1/1,verboseRMMH=T, maxRMMHiter=200)
#'
#' @export rmmh
#-------------------------------------------------------------------------------
#
# if(1==0){
#   chains=1
#   y = sdirt$y
#   objFun = dich_response_model
#   estOmega = T
#   estNu = T
#   estZeta = T
#   lambda = sdirt$lambda
#   gamma = sdirt$gamma
#   omega0 = array(data = 0, dim = dim(sdirt$omega))
#   nu0 = array(data = 0, dim = c(ncol(sdirt$nu), 1))
#   zeta0 = array(data = 0, dim = dim(sdirt$zeta))
#   omega_mu = sdirt$omega_mu
#   omega_sigma2 = sdirt$omega_sigma2
#   nu_mu = matrix(sdirt$nu_mu)
#   nu_sigma2 = matrix(sdirt$nu_sigma2)
#   zeta_mu = sdirt$zeta_mu
#   zeta_sigma2 = sdirt$zeta_sigma2
#   burn=50
#   thin=10
#   min_tune=10
#   tune_int=10
#   max_tune=100
#   niter=100
#   weight=1/1
#   verboseRMMH=T
#   maxRMMHiter=200
# }

rmmh <- function(
  chains = 1,
  y = y,
  objFun = NULL,
  estOmega = T,
  estNu = T,
  estZeta = T,
  lambda = lambda,
  gamma = gamma,
  omega0 = NULL,
  nu0 = NULL,
  zeta0 = NULL,
  omega_mu = NULL,
  omega_sigma2 = NULL,
  nu_mu = NULL,
  nu_sigma2 = NULL,
  zeta_mu = NULL,
  zeta_sigma2 = NULL,
  burn = NULL,
  thin = NULL,
  min_tune = NULL,
  tune_int = NULL,
  max_tune = NULL,
  niter = NULL,
  weight = NULL,
  verboseRMMH = T,
  maxRMMHiter = 200
) {
  # STEP 0: MCMC burnin
  if(verboseRMMH) {
    cat(
      "MCMC Burn-In Start Time",
      format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
      "\n",
      sep=" "
    )
  }
  mcmhburn=mcmhMC(
    chains=chains, y = y, objFun = objFun, estOmega = estOmega, estNu = estNu,
    estZeta = estZeta, lambda = lambda, gamma = gamma, omega0 = omega0,
    nu0 = nu0, zeta0 = zeta0, omega_mu = omega_mu, omega_sigma2 = omega_sigma2,
    nu_mu = nu_mu, nu_sigma2 = nu_sigma2, zeta_mu = zeta_mu,
    zeta_sigma2 = zeta_sigma2, burn = burn, thin = thin, min_tune = min_tune,
    tune_int = tune_int, max_tune = max_tune, niter = niter, weight = weight
  )

  # Update initial estimates and variance of candidates
  if(estOmega){
    omega0 <- mcmhburn$omegaEAP
#    omega_sigma2 <- mcmhburn$mcmhdraws$chain1$omega_sigma2
  }
  if(estNu){
    nu0 <- mcmhburn$nuEAP
#    nu_sigma2 <- mcmhburn$mcmhdraws$chain1$nu_sigma2
  }
  if(estZeta){
    zeta0 <- mcmhburn$zetaEAP
#    zeta_sigma2 <- mcmhburn$mcmhdraws$chain1$zeta_sigma2
  }

  # Set up RMMH loop parmaeters
  tol = .001
  ll0 <- Inf
  ll1 <- 0
  iter <- 0
  info0 <- diag(10, ncol(omega0))
  go <- T
  ongoingOmega <- omega0

  #start iteration
  if(verboseRMMH) {
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
  while(go && (iter < maxRMMHiter)) {
    iter <- iter + 1
    # STEP 1: Stochastic imputation
    mkDrawsAtIteratoinK = mcmhMC(
      chains = 3, y = y, objFun = objFun, estOmega = estOmega,
      estNu = estNu, estZeta = estZeta, lambda = lambda, gamma = gamma,
      omega0 = omega0,
      nu0 = nu0,
      zeta0 = zeta0,
      omega_mu = omega_mu, omega_sigma2 = omega_sigma2,
      nu_mu = matrix(nu_mu), nu_sigma2 = matrix(nu_sigma2),
      zeta_mu = zeta_mu, zeta_sigma2 = zeta_sigma2,
      burn=0, thin=10, min_tune=50, tune_int=50, max_tune = 1000, niter = 2000,
      weight=1/1
    )

    # STEP 2: Stochastic approximation
    # !!!in the orignal scrip this is a 1X16Xniter matrix
    grad <- NULL
    for(i in 1:(length(seq(from = burn + 1, to = niter, by = thin)) * chains)) {
      grad <- abind::abind(
        grad,
        #!!!in the orignal scrip this is a 1X16 matrix
        array(unlist(dich_response_deriv(
          y = y,
          nu = array(
            data = if(estNu) {
              array(data = abind::abind(lapply(mkDrawsAtIteratoinK$mcmhdraws,function(x){x[["nu_draws"]]}), along = 1)[i,,], dim = dim(nu0))
            } else {
              nu0
            },
            dim = dim(y)
          ),
          lambda = lambda,
          gamma = gamma,
          omega = omega0,
          zeta = if(estZeta) {
            matrix(data = abind::abind(lapply(mkDrawsAtIteratoinK$mcmhdraws,function(x){x[["zeta_draws"]]}), along = 1)[i,,], nrow = nrow(zeta0), ncol = ncol(zeta0), byrow = T)
          } else{
            zeta0
          },
          omega_mu = omega_mu,
          omega_sigma2 = omega_sigma2,
          zeta_mu = zeta_mu,
          zeta_sigma2 = zeta_sigma2
        )$fpd),c(nrow(omega0),ncol(omega0))),
        along=3
      )
    }
    # !!!in the orignal scrip this is a 1X16 matrix
    grad <- apply(X = grad, MARGIN = c(1,2), FUN = mean)
    # # !!!in the orignal scrip this is a 16X16 matrix
    info <- NULL
    for(i in 1:(length(seq(from = burn + 1, to = niter, by = thin)) * chains)) {
      info <- abind::abind(
        info,
        array(unlist(dich_response_deriv(
          y = y,
          nu = array(
            data = if(estNu) {
              array(data = abind::abind(lapply(mkDrawsAtIteratoinK$mcmhdraws,function(x){x[["nu_draws"]]}), along = 1)[i,,], dim = dim(nu0))
            } else {
              nu0
            },
            dim = dim(y)
          ),
          lambda = lambda,
          gamma = gamma,
          omega = omega0,
          zeta = if(estZeta) {
            matrix(data = abind::abind(lapply(mkDrawsAtIteratoinK$mcmhdraws,function(x){x[["zeta_draws"]]}), along = 1)[i,,], nrow = nrow(zeta0), ncol = ncol(zeta0), byrow = T)
          } else{
            zeta0
          },
          omega_mu = omega_mu,
          omega_sigma2 = omega_sigma2,
          zeta_mu = zeta_mu,
          zeta_sigma2 = zeta_sigma2
        )$post_info),c(ncol(omega0),ncol(omega0))),
        along=3
      )
    }
    info <- apply(X = info, MARGIN = c(1,2), FUN = mean)
    gain = 1 / iter
    info1 <- info0 + gain * (info - info0)
    # STEP 3: Robbins-Monro updat
    omega1 <- omega0 + gain * grad %*% solve(info1)
    p <-objFun(y = y, nu = array(data = nu0, dim = dim(y)), lambda = lambda,
                  gamma = gamma, omega = omega1, zeta = zeta0)$p
    ll1 <- sum(x = log((p^y) * (1 - p)^(1 - y)), na.rm = T)
    #Test for completion
    if(any(abs(omega1 - omega0) > tol) & iter < maxRMMHiter) {
      go = T
      #update values
      info0 <- info1
      omega0 <- omega1
      ongoingOmega <- c(ongoingOmega, omega1)
      if(estNu) {
        nu0 <- mkDrawsAtIteratoinK$nuEAP
      }
      if(estZeta) {
        zeta0 <- mkDrawsAtIteratoinK$zetaEAP
      }
      if(verboseRMMH) {
        cat(
          "\r",
          "... at iteration ",
          iter,
          " MAP omega is ",
          format(x = round(x = omega1, digits = 3), nsmall = 3),
          " logLik is ",
          format(x = round(x = ll1, digits = 4), nsmall  =4),
          sep = ""
          )
        }
    } else if(all(abs(omega1 - omega0) < tol) & iter < maxRMMHiter){
      go=F
      if(verboseRMMH) {
        cat(
          "\n... algorithm converged at",
          format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
          "\n",
          sep=" "
          )
        }
    } else if (all(abs(omega1-omega0)<tol) & iter == maxRMMHiter) {
      if(verboseRMMH) {
        cat(
          "\n... algorithm failed to converge at",
          format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
          "\n",
          sep=" "
          )
        }
    }
  }
  return(list("omega1" = omega1, "info1" = info1, "ll1" = ll1,
              "ongoingOmega" = ongoingOmega))
}
