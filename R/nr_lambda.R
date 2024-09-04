#-------------------------------------------------------------------------------
#' Newton Raphson Parameter Estimates for Lambda
#'
#' This function calculates Newton Raphson parameter estimates.
#'
#' @examples
#'
#' nr_lambda(y = ex1$y, omega0 = ex1$omega, gamma0 = ex1$gamma,
#'             lambda0 = ex1$lambda, zeta0 = ex1$zeta, nu0 = ex1$nu,
#'             kappa0 = ex1$kappa, lambda_mu = ex1$lambda_mu,
#'             lambda_sigma2 = ex1$lambda_sigma2, link  = "probit",
#'             obj_fun = dich_response_model)
#'
#' @keywords internal
#-------------------------------------------------------------------------------

nr_lambda <- function(y = NULL, omega0 = NULL, gamma0 = NULL, lambda0 = NULL,
                      zeta0 = NULL, nu0 = NULL, kappa0 = NULL, lambda_mu = NULL,
                      lambda_sigma2 = NULL, link  = "probit", obj_fun = NULL,
                      tol = 1e-9, max_iter = 1000, verbose = TRUE) {

  lambda1 <- lambda0
  kappa <- if (is.null(x = kappa)) {
    array(data = 0, dim = c(nrow(x = y), 1))
  } else {
    array(data = 1, dim = c(nrow(x = y), 1)) %*% t(kappa)
  }

  fx <- deriv_lambda(
    y = y,
    omega = omega0,
    gamma = gamma0,
    lambda = lambda1,
    zeta = zeta0,
    nu = nu0,
    kappa = kappa0,
    lambda_mu = lambda_mu,
    lambda_sigma2 = lambda_sigma2,
    link  = "probit"
  )

  iter <- 0

  while (any(abs(unlist(fx$fpd)) > tol) && (iter < max_iter)) {

    lambda1 <- t(t(sapply(
      X = 1:nrow(x = lambda1),
      FUN = function(x) {
        lambda1[x, ] - fx$fpd[[x]] %*% solve(fx$spd[[x]])
      }
    )))

    fx <- deriv_lambda(
      y = y,
      omega = omega0,
      gamma = gamma0,
      lambda = lambda1,
      zeta = zeta0,
      nu = nu0,
      kappa = kappa0,
      lambda_mu = lambda_mu,
      lambda_sigma2 = lambda_sigma2,
      link  = "probit"
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
  return("lambda1" = lambda1)
}
