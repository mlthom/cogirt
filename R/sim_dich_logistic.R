#' Simulate Data from a Dichotomous Logistic Response Model
#'
#' This function simulates data for a dichotmous logist response model. The
#' simulation example is currently designed to only simulate the SD-IRT model.
#'
#' @param N Sample size.
#' @param I Number of items per condition.
#' @param J Number of eperimental conditions.
#' @param M Number of ability (or trait) dimensions.
#' @param K Number of experimental contrasts (includes intercept).
#' @param nu_mu Mean of the item intercept parameters (scalar).
#' @param nu_sigma2 Variance of the item intercept parameters (scalar).
#' @param omega_mu Vector of means for the subject-level effects of the
#' experimental manipulation.
#' @param omega_sigma2 Covariance matrix for the subject-level effects of the
#' experimental manipulation.
#' @param zeta_mu Vector of means for the condition-level prediction errors.
#' @param zeta_sigma2 Covariance matrix for the condition-level prediction
#' errors.
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

sim_dich_logistic <- function(N, I, J, M, K, nu_mu, nu_sigma2, lambda, gamma, omega_mu, omega_sigma2, zeta_mu, zeta_sigma2, link  = 'logit') {
  nu <- matrix(data = rnorm(n = I*J, mean = nu_mu, sd = nu_sigma2), nrow = I * J, ncol = N)
  omega<- mvrnorm(n = N, mu = omega_mu, Sigma = omega_sigma2)
  zeta <- mvrnorm(n = N, mu = zeta_mu, Sigma = zeta_sigma2 * diag(x = 1, nrow = M * J))
  scaling_constant <- if(link == 'logit') {
    1.702
  } else if(link == 'probit'){
    1.000
  }
  epsilon <- matrix(data = rnorm(n = I * J * N, mean = 0, sd = scaling_constant), nrow = I * J, ncol = N)
  ystar <- t(nu + lambda %*% gamma %*% t(omega) + lambda %*% t(zeta) + epsilon)
  y <- apply(data.frame(ystar > 0), 2, as.numeric)
  colnames(y) <- c(sapply(X = 1:J, FUN = function(j){paste("item", 1:I, "cond", j, sep="")}))
  rownames(y) <- 1:N
  simdat <- list(y = y, nu = nu, lambda = lambda, gamma = gamma, omega = omega, zeta = zeta)
  return(simdat)
}

