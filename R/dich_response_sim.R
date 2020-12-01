#-------------------------------------------------------------------------------
#' Simulate Dichotomous Response Model
#'
#' This function simulates data for a dichotomous response model framed using
#' generalized latent variable modeling (GLVM; Skrondal & Rabe-Hesketh, 2004).
#' The example is based on the signal detection item response theory (IRT)
#' model (Thomas et al., 2018).
#'
#' @param I Number of items per condition.
#' @param J Number of conditions.
#' @param K Number of examinees
#' @param M Number of ability (or trait) dimensions.
#' @param N Number of contrasts (should include intercept).
#' @param nu_mu Mean of the item intercept parameters (scalar).
#' @param nu_sigma2 Variance of the item intercept parameters (scalar).
#' @param omega_mu Vector of means for the examinee-level effects of the
#' experimental manipulation (1 by MN).
#' @param omega_sigma2 Covariance matrix for the examinee-level effects of the
#' experimental manipulation (MN by MN).
#' @param zeta_mu Vector of means for the condition-level effects nested within
#' examinees (1 by JM).
#' @param zeta_sigma2 Covariance matrix for the condition-level effects nested
#' within examinees (JM by JM).
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
#' #Multiple Subjects
#'
#' I = 20
#' J = 5
#' K = 50
#' M = 2
#' N = 3
#' nu_mu = 0
#' nu_sigma2 = 0.2
#' omega_mu <- matrix(data = c(4.00, -2.00, 0.01, 0.00, 0.50, 0.00), nrow = 1,
#' ncol = M * N)
#' omega_sigma2 <- diag(x = c(0.50, 0.10, 0.01, 0.25, 0.10, 0.01), nrow = M * N)
#' zeta_mu <- matrix(data = rep(x = 0, times = M * J), nrow = 1, ncol = J * M)
#' zeta_sigma2 <- diag(x = 0.2, nrow = J * M, ncol = J * M)
#' item_type <- rbinom(n = I * J, size = 1, prob = .7) + 1
#' # Equation 12 Thomas et al. (2018)
#' measure_weights <-
#'   matrix(data = c(0.5, -1.0, 0.5, 1.0), nrow = 2, ncol = M, byrow = T)
#' lambda <- matrix(data = 0, nrow = I * J, ncol = J * M)
#' for(j in 1:J) {
#'   lambda[(1 + (j - 1) * I):(j * I), (1 + (j - 1) * M):(j * M)] <-
#'     measure_weights[item_type, ][(1 + (j - 1) * I):(j * I), ]
#' }
#' contrast_codes <- cbind(1, contr.poly(n = J))[, 1:N]
#' gamma <- matrix(data = 0, nrow = J * M, ncol = M * N)
#' for(j in 1:J) {
#'   for(m in 1:M) {
#'     gamma[(m + M * (j - 1)), (((m - 1) * N + 1):((m - 1) * N + N))] <-
#'     contrast_codes[j, ]
#'   }
#' }
#'
#'
#' sdirt <- dich_response_sim(I = I, J = J, K = K, M = M, N = N, nu_mu = nu_mu,
#'                   nu_sigma2 = nu_sigma2, lambda = lambda, gamma = gamma,
#'                   omega_mu = omega_mu, omega_sigma2 = omega_sigma2,
#'                   zeta_mu = zeta_mu, zeta_sigma2 = zeta_sigma2)
#'
#' #Single Subject
#'
#' I = 50
#' J = 5
#' K = 1
#' M = 2
#' N = 3
#' nu_mu = 0
#' nu_sigma2 = 0.2
#' omega_mu <- matrix(data = c(4.00, -2.00, 0.01, 0.00, 0.50, 0.00), nrow = 1,
#' ncol = M * N)
#' omega_sigma2 <- diag(x = c(0.50, 0.10, 0.01, 0.25, 0.10, 0.01), nrow = M * N)
#' zeta_mu <- matrix(data = rep(x = 0, times = M * J), nrow = 1, ncol = J * M)
#' zeta_sigma2 <- diag(x = 0.2, nrow = J * M, ncol = J * M)
#' item_type <- rbinom(n = I * J, size = 1, prob = .7) + 1
#' # Equation 12 Thomas et al. (2018)
#' measure_weights <-
#'   matrix(data = c(0.5, -1.0, 0.5, 1.0), nrow = 2, ncol = M, byrow = T)
#' lambda <- matrix(data = 0, nrow = I * J, ncol = J * M)
#' for(j in 1:J) {
#'   lambda[(1 + (j - 1) * I):(j * I), (1 + (j - 1) * M):(j * M)] <-
#'     measure_weights[item_type, ][(1 + (j - 1) * I):(j * I), ]
#' }
#' contrast_codes <- cbind(1, contr.poly(n = J))[, 1:N]
#' gamma <- matrix(data = 0, nrow = J * M, ncol = M * N)
#' for(j in 1:J) {
#'   for(m in 1:M) {
#'     gamma[(m + M * (j - 1)), (((m - 1) * N + 1):((m - 1) * N + N))] <-
#'     contrast_codes[j, ]
#'   }
#' }
#'
#'
#' sdirtSS <- dich_response_sim(I = I, J = J, K = K, M = M, N = N,
#'                   nu_mu = nu_mu, nu_sigma2 = nu_sigma2, lambda = lambda,
#'                   gamma = gamma, omega_mu = omega_mu,
#'                   omega_sigma2 = omega_sigma2, zeta_mu = zeta_mu,
#'                   zeta_sigma2 = zeta_sigma2)
#'
#' @export dich_response_sim
#-------------------------------------------------------------------------------

dich_response_sim <- function(I, J, K,  M, N, nu_mu, nu_sigma2, lambda, gamma,
                              omega_mu, omega_sigma2, zeta_mu, zeta_sigma2,
                              link  = "logit") {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package \"MASS\" needed for the dich_response_sim function to work.
         Please install.",
         call. = FALSE)
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package \"MASS\" needed for the dich_response_sim function to work.
         Please install.",
         call. = FALSE)
  }
  nu <- matrix(
    data = MASS::mvrnorm(
      n = I * J,
      mu = nu_mu,
      Sigma = nu_sigma2
    ),
    nrow = K,
    ncol = I * J,
    byrow = T
  )
  omega <-  matrix(
    data = MASS::mvrnorm(n = K, mu = omega_mu, Sigma = omega_sigma2),
    nrow = K,
    ncol = N * M,
    byrow = F
  )
  zeta <- matrix(
    data = MASS::mvrnorm(
      n = K,
      mu = zeta_mu,
      Sigma = zeta_sigma2 * diag(x = 1, nrow = M * J)
    ),
    nrow = K,
    ncol = J * M,
    byrow = F
  )
  # 'scaling_constant' used to define epsilon variance
  scaling_constant <- if (link == "logit") {
    1.702
  } else if (link == "probit") {
    1.000
  }
  epsilon <- matrix(
    data = MASS::mvrnorm(n = I * J * K, mu = 0, Sigma = scaling_constant),
    nrow = K,
    ncol = I * J)
  ystar <- nu + omega %*% t(gamma) %*% t(lambda) + zeta %*% t(lambda) + epsilon
  y <- apply(data.matrix(ystar > 0), 1:2, as.numeric)
  colnames(y) <- c(
    sapply(
      X = 1:J,
      FUN = function(j) {
        paste("item", 1:I, "cond", j, sep = "")
      }
    )
  )
  rownames(y) <- paste("examinee", 1:K, sep = "")
  condition <- c(sapply(X = 1:J, FUN = rep, I * J / J))
  simdat <- list("y" = y, "ystar" = ystar, "nu" = nu, "lambda" = lambda,
                 "gamma" = gamma, "omega" = omega, "zeta" = zeta,
                 "nu_mu" = nu_mu, "nu_sigma2" = nu_sigma2,
                 "omega_mu" = omega_mu,
                 "omega_sigma2" = omega_sigma2, "zeta_mu" = zeta_mu,
                 "zeta_sigma2" = zeta_sigma2, "condition" = condition,
                 "item_type" = item_type
  )
  return(simdat)
}
