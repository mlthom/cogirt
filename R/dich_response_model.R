#-------------------------------------------------------------------------------
#' Dichotomous Response Model
#'
#' This function calculates predictions and log-likelihood values for a
#' dichotomous response model framed using generalized latent variable modeling
#' (GLVM; Skrondal & Rabe-Hesketh, 2004).
#'
#' @param y Matrix of item responses (N by I * J).
#' @param nu Matrix of item intercept parameters (N by I * J).
#' @param lambda Matrix of item structure parameters (I * J by J * M).
#' @param gamma Matrix of experimental structure parameters (J * M by K * M).
#' @param omega Examinee-level effects of the experimental manipulation
#' (N by K * M).
#' @param zeta Condition-level effects nested within examineess (N by J * M).
#' @param link Choose between logit or probit link functions.
#'
#' @return p = response probability matrix (N by I * J); yhatstar = latent
#' response variate matrix (N by I * J); loglikelihood = model log-likelihood
#' (scalar).
#'
#' @section Dimensions:
#' N = Number of examinees; I = Number of items per condition; J = Number of
#' conditions; M = Number of ability (or trait) dimensions; K = Number of
#' contrasts.
#'
#' @references
#'
#' Skrondal, A., & Rabe-Hesketh, S. (2004). \emph{Generalized latent variable
#' modeling: Multilevel, longitudinal, and structural equation models}. Boca
#' Raton: Chapman & Hall/CRC.
#'
#' @examples
#' mod <- dich_response_model(y = sdirt$y, nu = sdirt$nu, lambda = sdirt$lambda,
#'                      gamma = sdirt$gamma, omega = sdirt$omega, sdirt$zeta)
#'
#' @export dich_response_model
#-------------------------------------------------------------------------------

dich_response_model <- function(y, nu, lambda, gamma, omega, zeta, link  = 'logit') {
  yhatstar <- nu + omega %*% t(gamma) %*% t(lambda) + zeta %*% t(lambda)
  p <- if(link == 'logit') {
    plogis(yhatstar)
  } else if(link == 'probit'){
    pnorm(yhatstar)
  }
  ll <- sum(log((p^y) * (1 - p)^(1 - y)))
  mod <- list(p = p, yhatstar = yhatstar, loglikelihood = ll)
  return(mod)
}
