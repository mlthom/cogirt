#-------------------------------------------------------------------------------
#' Dichotomous Response Model
#'
#' This function calculates predictions and log-likelihood values for a
#' dichotomous response model framed using generalized latent variable modeling
#' (GLVM; Skrondal & Rabe-Hesketh, 2004).
#'
#' @param y Matrix of item responses (K by IJ).
#' @param omega Examinee-level effects of the experimental manipulation (K by
#' MN).
#' @param gamma Matrix of experimental structure parameters (JM by MN).
#' @param lambda Matrix of item structure parameters (IJ by JM).
#' @param zeta Condition-level effects of the experimental manipulation (K by
#' JM).
#' @param nu Matrix of item intercept parameters (IJ by 1).
#' @param kappa Matrix of item guessing parameters (IJ by 1).
#' @param link Choose between "logit" or "probit" link functions.
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
#' @export dich_response_model
#-------------------------------------------------------------------------------

dich_response_model <- function(y = NULL, omega = NULL, gamma = NULL,
                                lambda = NULL, zeta = NULL, nu = NULL,
                                kappa = NULL, link  = NULL) {
  link <- if (is.null(x = link)) {
    "probit"
  } else {
    link
  }
  kappa_mat <- if (is.null(x = kappa)) {
    array(data = 0, dim = dim(x = y))
  } else {
    array(data = 1, dim = c(nrow(x = y), 1)) %*% t(kappa)
  }
  nu_mat <- array(data = 1, dim = c(nrow(x = y), 1)) %*% t(nu)
  yhatstar <- nu_mat + omega %*% t(gamma) %*% t(lambda) + zeta %*% t(lambda)
  p <- if (link == "logit") {
    kappa_mat + (1 - kappa_mat) * plogis(q = yhatstar)
  } else if (link == "probit") {
    kappa_mat + (1 - kappa_mat) * pnorm(q = yhatstar)
  }
  ll <- sum(log(x = (p ^ y) * (1 - p) ^ (1 - y)))
  mod <- list(p = p, yhatstar = yhatstar, loglikelihood = ll)
  return(mod)
}
