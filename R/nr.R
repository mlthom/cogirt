#-------------------------------------------------------------------------------
#' Newton Raphson Parameter Estimates
#'
#' This function calculates Newton Raphson parameter estimates. It is used
#' mostly for testing.
#'
#' @export nr
#-------------------------------------------------------------------------------

nr <- function(rda = NULL, tol=1e-9, max_iter=100, verbose=T) {

  omegaNR <- rda$omega
  fx <- dich_response_deriv(
    y = rda$y,
    nu = rda$nu,
    lambda = rda$lambda,
    gamma = rda$gamma,
    omega = rda$omega,
    zeta = rda$zeta,
    omega_mu = rda$omega_mu,
    omega_sigma2 = rda$omega_sigma2,
    zeta_mu = rda$zeta_mu,
    zeta_sigma2 = rda$zeta_sigma2
  )
  iter <- 0
  while(any(abs(fx$fpd[[1]]) > tol) && (iter < max_iter)){

    omegaNR <- omegaNR - fx$fpd[[1]] %*% solve(fx$spd[[1]])

    #bounded solution to help convergence
    #if(any(abs(omegaNR)>(diag(rda$omega_sigma2)*5))){
    #  omegaNR[which(abs(omegaNR)>(diag(rda$omega_sigma2)*5))] <- (rnorm(length(omegaNR)) + sign(omegaNR)*(diag(rda$omega_sigma2)*5))[which(abs(omegaNR)>(diag(rda$omega_sigma2)*5))]
    #}

    fx <- dich_response_deriv(
      y = rda$y,
      nu = rda$nu,
      lambda = rda$lambda,
      gamma = rda$gamma,
      omega = omegaNR,
      zeta = rda$zeta,
      omega_mu = rda$omega_mu,
      omega_sigma2 = rda$omega_sigma2,
      zeta_mu = rda$zeta_mu,
      zeta_sigma2 = rda$zeta_sigma2
    )
    iter <- iter + 1
  }
  if(verbose==T){
    if(any(abs(fx$fpd[[1]]) > tol)){
      cat("Algorithm failed to converge\n")
    } else {
      cat("Algorithm converged\n")
    }
  }
  return("omegaNR" = omegaNR)
}

