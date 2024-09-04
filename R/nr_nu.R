#-------------------------------------------------------------------------------
#' Newton Raphson Parameter Estimates for Nu
#'
#' This function calculates Newton Raphson parameter estimates.
#'
#' @examples
#'
#' nr_nu(y = ex1$y, omega0 = ex1$omega, gamma0 = ex1$gamma,
#'             lambda0 = ex1$lambda, zeta0 = ex1$zeta, nu0 = ex1$nu,
#'             kappa0 = ex1$kappa, nu_mu = ex1$nu_mu,
#'             nu_sigma2 = ex1$nu_sigma2, link  = "probit",
#'             obj_fun = dich_response_model)
#'
#'@keywords internal
#-------------------------------------------------------------------------------

nr_nu <- function(y = NULL, omega0 = NULL, gamma0 = NULL, lambda0 = NULL,
                  zeta0 = NULL, nu0 = NULL, kappa0 = NULL, nu_mu = NULL,
                  nu_sigma2 = NULL, link  = "probit", obj_fun = NULL,
                  tol = 1e-9, max_iter = 1000, verbose = TRUE) {

  nu1 <- nu0
  kappa0 <- if (is.null(x = kappa)) {
    array(data = 0, dim = c(nrow(x = y), 1))
  } else {
    array(data = 1, dim = c(nrow(x = y), 1)) %*% t(kappa)
  }

  fx <- deriv_nu(
    y = y,
    omega = omega0,
    gamma = gamma0,
    lambda = lambda0,
    zeta = zeta0,
    nu = nu1,
    kappa = kappa0,
    nu_mu = nu_mu,
    nu_sigma2 = nu_sigma2
  )

  iter <- 0

  while (any(abs(unlist(fx$fpd)) > tol) && (iter < max_iter)) {

    nu1 <- t(t(sapply(
      X = 1:nrow(x = nu1),
      FUN = function(x) {
        nu1[x, 1] - fx$fpd[[x]] %*% solve(fx$spd[[x]])
      }
    )))

    nu1 <- nu1

    fx <- deriv_nu(
      y = y,
      omega = omega0,
      gamma = gamma0,
      lambda = lambda0,
      zeta = zeta0,
      nu = nu1,
      kappa = kappa0,
      nu_mu = nu_mu,
      nu_sigma2 = nu_sigma2
    )
    iter <- iter + 1
  }
  if (verbose == TRUE) {
    if (any(abs(unlist(fx$fpd)) > tol)) {
      cat("Algorithm failed to converge\n")
    } else {
      cat("Algorithm converged\n")
    }
  }
  return("nu1" = nu1)
}
