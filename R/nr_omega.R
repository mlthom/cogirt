#-------------------------------------------------------------------------------
#' Newton Raphson Parameter Estimates for Omega Parameters
#'
#' This function calculates Newton Raphson parameter estimates.
#'
#' @examples
#'
#' nr_omega(y = ex1$y, omega0 = ex1$omega, gamma0 = ex1$gamma,
#'             lambda0 = ex1$lambda, zeta0 = ex1$zeta, nu0 = ex1$nu,
#'             kappa0 = ex1$kappa, omega_mu = ex1$omega_mu,
#'             omega_sigma2 = ex1$omega_sigma2, zeta_mu = ex1$zeta_mu,
#'             zeta_sigma2 = ex1$zeta_sigma2, link = "probit",
#'             obj_fun = dich_response_model)
#'
#' @keywords internal
#-------------------------------------------------------------------------------

nr_omega <- function(y = NULL, omega0 = NULL, gamma0 = NULL, lambda0 = NULL,
                     zeta0 = NULL, nu0 = NULL, kappa0 = NULL, omega_mu = NULL,
                     omega_sigma2 = NULL, zeta_mu = NULL,
                     zeta_sigma2 = NULL, link  = "probit", obj_fun = NULL,
                     tol = 1e-9, max_iter = 1000, verbose = TRUE) {

  omega1 <- omega0

  kappa0 <- if (is.null(x = kappa)) {
    array(data = 0, dim = c(nrow(x = y), 1))
  } else {
    array(data = 1, dim = c(nrow(x = y), 1)) %*% t(kappa)
  }

  fx <- deriv_omega(
    y = y,
    omega = omega1,
    gamma = gamma0,
    lambda = lambda0,
    zeta = zeta0,
    nu = nu0,
    kappa = kappa0,
    omega_mu = omega_mu,
    omega_sigma2 = omega_sigma2,
    zeta_mu = zeta_mu,
    zeta_sigma2 = zeta_sigma2
  )

  iter <- 0

  while (any(abs(unlist(fx$fpd)) > tol) && (iter < max_iter)) {

    omega1 <- t(t(sapply(
      X = 1:nrow(x = omega1),
      FUN = function(x) {
        omega1[x, ] - fx$fpd[[x]] %*% solve(fx$spd[[x]])
      }
    )))

    fx <- deriv_omega(
      y = y,
      omega = omega1,
      gamma = gamma0,
      lambda = lambda0,
      zeta = zeta0,
      nu = nu0,
      kappa = kappa0,
      omega_mu = omega_mu,
      omega_sigma2 = omega_sigma2,
      zeta_mu = zeta_mu,
      zeta_sigma2 = zeta_sigma2
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
  return("omega1" = omega1)
}
