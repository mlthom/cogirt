#-------------------------------------------------------------------------------
#' Dichotomous Response Model
#'
#' This function calculates predictions and log-likelihood values for a
#' dichotomous response model framed using generalized latent variable modeling
#' (GLVM; Skrondal & Rabe-Hesketh, 2004).
#'
#' @param y Matrix of item responses (K by IJ).
#' @param nu Matrix of item intercept parameters (K by IJ).
#' @param omega Examinee-level effects of the experimental manipulation (K by
#' M*N).
#' @param gamma Matrix of experimental structure parameters (JM by MN).
#' @param lambda Matrix of item structure parameters (IJ by JM).
#' @param zeta Condition-level effects of the experimental manipulation (K by
#' JM).
#' @param link Choose between logit or probit link functions.
#'
#' @return p = response probability matrix (K by IJ); yhatstar = latent response
#' variate matrix (K by IJ); loglikelihood = model log-likelihood (scalar).
#'
#' @section Dimensions:
#' I = Number of items per condition; J = Number of conditions; K = Number of
#' examinees; M Number of ability (or trait) dimensions; N Number of contrasts
#' (should include intercept).
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

dich_response_model <- function(y, nu, lambda, gamma, omega, zeta,
                                link  = 'logit') {
  yhatstar <- nu + omega %*% t(gamma) %*% t(lambda) + zeta %*% t(lambda)
  p <- if(link == 'logit') {
    plogis(q = yhatstar)
  } else if(link == 'probit'){
    pnorm(q = yhatstar)
  }
  ll <- sum(log((p^y) * (1 - p)^(1 - y)))
  mod <- list(p = p, yhatstar = yhatstar, loglikelihood = ll)
  return(mod)
}
