#-------------------------------------------------------------------------------
#' Simulate Dichotomous Response Model
#'
#' This function simulates data for a dichotomous response model framed using
#' generalized latent variable modeling (GLVM; Skrondal & Rabe-Hesketh, 2004).
#' The example is based on the signal detection item response theory (IRT)
#' model (Thomas et al., 2018).
#'
#' @param N Number of examinees
#' @param I Number of items per condition.
#' @param J Number of conditions.
#' @param M Number of ability (or trait) dimensions.
#' @param K Number of contrasts (should include intercept).
#' @param nu_mu Mean of the item intercept parameters (scalar).
#' @param nu_sigma2 Variance of the item intercept parameters (scalar).
#' @param omega_mu Vector of means for the examinee-level effects of the
#' experimental manipulation (1 by K * M).
#' @param omega_sigma2 Covariance matrix for the examinee-level effects of the
#' experimental manipulation (K * M by K * M).
#' @param zeta_mu Vector of means for the condition-level effects nested within
#' examinees (1 by J * M).
#' @param zeta_sigma2 Covariance matrix for the condition-level effects nested
#' within examinees (J * M by J * M).
#' @param link Choose between logit or probit links (i.e., sets epsilon variance
#' to 1.702 or 1.000 respectively).
#'
#' @return y = simulated response matrix; yhatstar = simulated latent response
#' probability matrix; [simulation_parameters]
#'
#' @references
#'
#' Skrondal, A., & Rabe-Hesketh, S. (2004). \emph{Generalized latent variable
#' modeling: Multilevel, longitudinal, and structural equation models}. Boca
#' Raton: Chapman & Hall/CRC.
#'
#' Thomas, M. L., Brown, G. G., Gur, R. C., Moore, T. M., Patt,
#' V. M., Risbrough, V. B., & Baker, D. G. (2018). A signal detection-item
#' response theory model for evaluating neuropsychological measures.
#' \emph{Journal of Clinical and Experimental Neuropsychology, 40(8)}, 745-760.
#'
#' @examples
#' N = 50
#' I = 200
#' J = 5
#' M = 2
#' K = 2
#' nu_mu = 0
#' nu_sigma2 = 0.2
#' omega_mu <- matrix(data = c(4, -0.5, 0, .5), nrow = 1, ncol = K * M)
#' omega_sigma2 <- diag(x = c(.5, 0.1, 0.25, 0.1), nrow = M * K)
#' zeta_mu <- matrix(data = rep(x = 0, times = M * J), nrow = 1, ncol = J * M)
#' zeta_sigma2 <- diag(x = 0.2, nrow = J * M, ncol = J * M)
#' item_type <- rbinom(n = I * J, size = 1, prob = .7) + 1
#' measure_weights <-
#'   matrix(data = c(0.5, -1.0, 0.5, 1.0), nrow = 2, ncol = 2, byrow = T)
#' lambda <- matrix(data = 0, nrow =I * J, ncol = J * M)
#' for(j in 1:J){
#'   lambda[(1 + (j - 1) * I):(j * I), (1 + (j - 1) * M):(j * M)] <-
#'     measure_weights[item_type, ][(1 + (j - 1) * I):(j * I), ]
#' }
#'
#' contrast_codes <- cbind(1, contr.poly(n = J))[, 1:K]
#' gamma <- matrix(data = 0, nrow = J * M, ncol = K * M)
#' for(j in 1:J){
#'   for(m in 1:M){
#'     gamma[m + M * (j - 1), (1 + (m - 1) * M):(m * M)] <- contrast_codes[j, ]
#'   }
#' }
#'
#'
#' sdirt <- dich_response_sim(N = N, I = I, J = J, M = M, K = K, nu_mu = nu_mu,
#'                   nu_sigma2 = nu_sigma2, lambda = lambda, gamma = gamma,
#'                   omega_mu = omega_mu, omega_sigma2 = omega_sigma2,
#'                   zeta_mu = zeta_mu, zeta_sigma2 = zeta_sigma2)
#' @export dich_response_sim
#-------------------------------------------------------------------------------

dich_response_sim <- function(N, I, J, M, K, nu_mu, nu_sigma2, lambda, gamma,
                              omega_mu, omega_sigma2, zeta_mu, zeta_sigma2,
                              link  = 'logit') {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package \"MASS\" needed for this function to work. Please install.",
         call. = FALSE)
  }
  # Matrices used to prevent dimension dropping
  nu <- matrix(data = MASS::mvrnorm(n = I * J,
                                    mu = nu_mu,
                                    Sigma = nu_sigma2),
               nrow = N,
               ncol = I * J,
               byrow = T)
  omega <-  matrix(data = MASS::mvrnorm(n = N,
                                        mu = omega_mu,
                                        Sigma = omega_sigma2),
                   nrow = N,
                   ncol = K * M,
                   byrow = F)
  zeta <- matrix(data = MASS::mvrnorm(n = N,
                                      mu = zeta_mu,
                                      Sigma = zeta_sigma2 *
                                        diag(x = 1, nrow = M * J)),
                 nrow = N,
                 ncol = J * M,
                 byrow = F)
  # 'scaling_constant' used to define epsilon variance
  scaling_constant <- if(link == 'logit') {
    1.702
  } else if(link == 'probit'){
    1.000
  }
  epsilon <- matrix(data = MASS::mvrnorm(n = I * J * N,
                                         mu = 0,
                                         Sigma = scaling_constant),
                    nrow = N,
                    ncol = I * J)
  ystar <- nu + omega %*% t(gamma) %*% t(lambda) + zeta %*% t(lambda) + epsilon
  y <- apply(data.matrix(ystar > 0), 1:2, as.numeric)
  colnames(y) <- c(sapply(X = 1:J,
                          FUN = function(j){
                            paste("item", 1:I, "cond", j, sep="")
                          }))
  rownames(y) <- paste("examinee", 1:N, sep = "")
  simdat <- list(y = y, ystar = ystar, nu = nu, lambda = lambda, gamma = gamma,
                 omega = omega, zeta = zeta, omega_mu = omega_mu,
                 omega_sigma2 = omega_sigma2, zeta_mu = zeta_mu,
                 zeta_sigma2 = zeta_sigma2)
  return(simdat)
}
