#' Dichotomous Response Model
#'
#' This function produces predictions and log-likelihood values for a dichotmous
#' response model.
#'
#' @param y Matrix of item responses (I * J by N).
#' @param nu Matrix of item intercept parameters (I * J by N).
#' @param lambda Matrix of item structure parameters (I * J by J * M).
#' @param gamma Matrix of experimental structure parameters (J * M by K * M).
#' @param omega Subject-level effects of the experimental manipulation (K * M by N).
#' @param zeta Condition-level prediction errors (N by J * M).
#' @param link Choose between logit or probit link functions.
#' @return yhat = response probability; yhatstar = latent response probability;
#' loglikelihood = model log-likelihood.
#' @examples
#' mod <- dich_response(y = sdirt$y, nu = sdirt$nu, lambda = sdirt$lambda,
#'                      gamma = sdirt$gamma, omega = sdirt$omega, sdirt$zeta)
#' @export

dich_response <- function(y, nu, lambda, gamma, omega, zeta, link  = 'logit') {
  yhatstar <- t(nu + lambda %*% gamma %*% t(omega) + lambda %*% t(zeta))
  yhat <- if(link == 'logit') {
    plogis(yhatstar)
  } else if(link == 'probit'){
    pnorm(yhatstar)
  }
  ll <- sum(log((yhat^y) * (1 - yhat)^(1 - y)))
  mod <- list(yhat = yhat, yhatstar = yhatstar, loglikelihood = ll)
  return(mod)
}
