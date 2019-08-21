#' Dichotomous Logistic Response Model
#'
#' This function produces predictions and log-likelihood values for a dichotmous
#' logist response model.
#'
#' @param y Vector of item responses.
#' @param nu Vector of item intercept parameters.
#' @param lambda Matrix of item structure parameters.
#' @param gamma Matrix of experimental structure parameters.
#' @param omega Subject-level effects of the experimental manipulation.
#' @param zeta Condition-level prediction errors.
#' @param link choose between logit and probit models.
#' @examples
#' N = 50
#' I = 20
#' J = 5
#' M = 2
#' K = 2
#' nu_mu = 0
#' nu_sigma2 = 0.2
#' zeta_mu <- rep(x = 0, times = M * J)
#' zeta_sigma2 <- 0.2
#' omega_mu <- c(4, -0.5, 0, .5)
#' omega_sigma2 <- diag(x = c(.5, 0.1, 0.25, 0.1), nrow = M * K)
#' item_type <- rbinom(n = I*J, size = 1, prob = .7) + 1
#' measure_weights <-
#'   matrix(data = c(0.5, -1.0, 0.5, 1.0), nrow = 2, ncol = 2, byrow = T)
#' lambda <- matrix(data = 0, nrow =I * J, ncol =M * J)
#' for(j in 1:J){
#'   lambda[(1 + (j - 1) * I):(j * I), (1 + (j - 1) * M):(j * M)] <-
#'     measure_weights[item_type, ][(1 + (j - 1) * I):(j * I), ]
#' }
#'
#' contrast_codes <- cbind(1, contr.poly(n = J))[, 1:K]
#' gamma <- matrix(data = 0, nrow = M * J, ncol = K * M)
#' for(j in 1:J){
#'   for(m in 1:M){
#'     gamma[m + M * (j - 1), (1 + (m - 1) * M):(m * M)] <- contrast_codes[j, ]
#'   }
#' }
#'
#'
#' sim_dich_logistic(N = N, I = I, J = J, M = M, K = K, nu_mu = nu_mu,
#'                   nu_sigma2 = nu_sigma2, lambda = lambda, gamma = gamma,
#'                   omega_mu = omega_mu, omega_sigma2 = omega_sigma2,
#'                   zeta_mu = zeta_mu, zeta_sigma2 = zeta_sigma2)
#' mod <- dich_logistic(y = simdat$y, nu = simdat$nu, lambda = simdat$lambda,
#'                      gamma = simdat$gamma, omega = simdat$omega, simdat$zeta)

dich_logistic <- function(y, nu, lambda, gamma, omega, zeta, link  = 'logit') {
  yhatstar <- t(nu + lambda %*% gamma %*% t(omega) + lambda %*% t(zeta))
  yhat <- if(link == 'logit') {
    plogis(yhatstar)
  } else if(link == 'probit'){
    pnorm(yhatstar)
  }
  ll <- sum(log((yhat^y) * (1 - yhat)^(1 - y)))
  mod <- list(predictions = yhat, loglikelihood = ll)
  return(mod)
}
