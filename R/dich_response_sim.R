#-------------------------------------------------------------------------------
#' Simulate Dichotomous Response Model
#'
#' This function calculates the matrix of first partial derivatives, the matrix
#' of second partial derivatives, and the information matrix for the posterior
#' distribution with respect to theta (ability) based on theslope-intercept form
#' of the item response theory model.
#'
#' @param I Number of items per condition.
#' @param J Number of conditions.
#' @param K Number of examinees
#' @param M Number of ability (or trait) dimensions.
#' @param N Number of contrasts (should include intercept).
#' @param omega Examinee-level effects of the experimental manipulation (K by
#' MN).
#' @param omega_mu Vector of means for the examinee-level effects of the
#' experimental manipulation (1 by MN).
#' @param omega_sigma2 Covariance matrix for the examinee-level effects of the
#' experimental manipulation (MN by MN).
#' @param gamma Matrix of experimental structure parameters (JM by MN).
#' @param lambda Matrix of item slope parameters (IJ by JM).
#' @param lambda_mu Vector of means for the item slope parameters (1 by JM)
#' @param lambda_sigma2 Covariance matrix for the item slope parameters (JM
#' by JM)
#' @param nu Matrix of item intercept parameters (K by IJ).
#' @param nu_mu Mean of the item intercept parameters (scalar).
#' @param nu_sigma2 Variance of the item intercept parameters (scalar).
#' @param zeta Condition-level effects of the experimental manipulation (K by
#' JM).
#' @param zeta_mu Vector of means for the condition-level effects nested within
#' examinees (1 by JM).
#' @param zeta_sigma2 Covariance matrix for the condition-level effects nested
#' within examinees (JM by JM).
#' @param kappa kappa	Matrix of item guessing parameters (IJ by 1). If kappa is not
#' provided, parameter values are set to 0.
#' @param link Choose between logit or probit link functions.
#' @param key Option key where  1 indicates target and 2 indicates distractor.
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
#'
#' # Example 1
#'
#' I <- 100
#' J <- 1
#' K <- 250
#' M <- 1
#' N <- 1
#' omega_mu <- matrix(data = 0, nrow = 1, ncol = M * N)
#' omega_sigma2 <- diag(x = 1, nrow = M * N)
#' gamma <- diag(x = 1, nrow = J * M, ncol = M * N)
#' lambda_mu <- matrix(data = 1, nrow = 1, ncol = M)
#' lambda_sigma2 <- diag(x = 0.25, nrow = M)
#' zeta_mu <- matrix(data = rep(x = 0, times = M * J), nrow = 1, ncol = J * M)
#' zeta_sigma2 <- diag(x = 0, nrow = J * M, ncol = J * M)
#' nu_mu <- matrix(data = 0, nrow = 1, ncol = 1)
#' nu_sigma2 <- matrix(data = 1, nrow = 1, ncol = 1)
#' set.seed(624)
#' ex1 <- dich_response_sim(I = I, J = J, K = K, M = M, N = N,
#'                          omega_mu = omega_mu, omega_sigma2 = omega_sigma2,
#'                          gamma = gamma, lambda_mu = lambda_mu,
#'                          lambda_sigma2 = lambda_sigma2, nu_mu = nu_mu,
#'                          nu_sigma2 = nu_sigma2, zeta_mu = zeta_mu,
#'                          zeta_sigma2 = zeta_sigma2)
#'
#' # Example 2
#'
#' I <- 100
#' J <- 1
#' K <- 50
#' M <- 2
#' N <- 1
#' omega_mu <- matrix(data = c(3.50, 1.00), nrow = 1, ncol = M * N)
#' omega_sigma2 <- diag(x = c(0.90, 0.30), nrow = M * N)
#' gamma <- diag(x = 1, nrow = J * M, ncol = M * N)
#' key <- rbinom(n = I * J, size = 1, prob = .7) + 1
#' measure_weights <-
#'   matrix(data = c(0.5, -1.0, 0.5, 1.0), nrow = 2, ncol = M, byrow = TRUE)
#' lambda <- matrix(data = 0, nrow = I * J, ncol = J * M)
#' for(j in 1:J) {
#'   lambda[(1 + (j - 1) * I):(j * I), (1 + (j - 1) * M):(j * M)] <-
#'     measure_weights[key, ][(1 + (j - 1) * I):(j * I), ]
#' }
#' zeta_mu <- matrix(data = rep(x = 0, times = M * J), nrow = 1, ncol = J * M)
#' zeta_sigma2 <- diag(x = 0, nrow = J * M, ncol = J * M)
#' nu_mu <- matrix(data = 0, nrow = 1, ncol = 1)
#' nu_sigma2 <- matrix(data = .2, nrow = 1, ncol = 1)
#' set.seed(624)
#' ex2 <- dich_response_sim(I = I, J = J, K = K, M = M, N = N,
#'                          omega_mu = omega_mu, omega_sigma2 = omega_sigma2,
#'                          gamma = gamma, lambda = lambda, nu_mu = nu_mu,
#'                          nu_sigma2 = nu_sigma2, zeta_mu = zeta_mu,
#'                          zeta_sigma2 = zeta_sigma2, key = key)
#'
#' # Example 3
#'
#' I <- 20
#' J <- 10
#' K <- 50
#' M <- 2
#' N <- 2
#' omega_mu <- matrix(data = c(2.50, -2.00, 0.50, 0.00), nrow = 1, ncol = M * N)
#' omega_sigma2 <- diag(x = c(0.90, 0.70, 0.30, 0.10), nrow = M * N)
#' contrast_codes <- cbind(1, contr.poly(n = J))[, 1:N]
#' gamma <- matrix(data = 0, nrow = J * M, ncol = M * N)
#' for(j in 1:J) {
#'   for(m in 1:M) {
#'     gamma[(m + M * (j - 1)), (((m - 1) * N + 1):((m - 1) * N + N))] <-
#'       contrast_codes[j, ]
#'   }
#' }
#' key <- rbinom(n = I * J, size = 1, prob = .7) + 1
#' measure_weights <-
#'   matrix(data = c(0.5, -1.0, 0.5, 1.0), nrow = 2, ncol = M, byrow = TRUE)
#' lambda <- matrix(data = 0, nrow = I * J, ncol = J * M)
#' for(j in 1:J) {
#'   lambda[(1 + (j - 1) * I):(j * I), (1 + (j - 1) * M):(j * M)] <-
#'     measure_weights[key, ][(1 + (j - 1) * I):(j * I), ]
#' }
#' zeta_mu <- matrix(data = rep(x = 0, times = M * J), nrow = 1, ncol = J * M)
#' zeta_sigma2 <- diag(x = 0.2, nrow = J * M, ncol = J * M)
#' nu_mu <- matrix(data = c(0.00), nrow = 1, ncol = 1)
#' nu_sigma2 <- matrix(data = c(0.20), nrow = 1, ncol = 1)
#' set.seed(624)
#' ex3 <- dich_response_sim(I = I, J = J, K = K, M = M, N = N,
#'                          omega_mu = omega_mu, omega_sigma2 = omega_sigma2,
#'                          gamma = gamma, lambda = lambda, nu_mu = nu_mu,
#'                          nu_sigma2 = nu_sigma2, zeta_mu = zeta_mu,
#'                          zeta_sigma2 = zeta_sigma2, key = key)
#'
#' # Example 4
#'
#' I <- 30
#' J <- 2
#' K <- 500
#' M <- 1
#' N <- 2
#' omega_mu <- matrix(data = c(0, -1), nrow = 1, ncol = M * N)
#' omega_sigma2 <- diag(x = c(1.00, 0.25), nrow = M * N)
#' contrast_codes <- cbind(1, contr.treatment(n = J))[, 1:N]
#' gamma <- matrix(data = 0, nrow = J * M, ncol = M * N)
#' for(j in 1:J) {
#'   for(m in 1:M) {
#'     gamma[(m + M * (j - 1)), (((m - 1) * N + 1):((m - 1) * N + N))] <-
#'     contrast_codes[j, ]
#'   }
#' }
#' lambda <- matrix(data = 0, nrow = I * J, ncol = J * M)
#' lam_vals <- rnorm(I, 1.5, .23)
#' for (j in 1:J) {
#'   lambda[(1 + (j - 1) * I):(j * I), (1 + (j - 1) * M):(j * M)] <- lam_vals
#' }
#' zeta_mu <- matrix(data = rep(x = 0, times = M * J), nrow = 1, ncol = J * M)
#' zeta_sigma2 <- diag(x = 0.2, nrow = J * M, ncol = J * M)
#' nu <- matrix(data = rnorm(n = I, mean = 0, sd = .25), nrow = I * J, ncol = 1)
#' set.seed(624)
#' ex4 <- dich_response_sim(I = I, J = J, K = K, M = M, N = N,
#'                          omega_mu = omega_mu, omega_sigma2 = omega_sigma2,
#'                          gamma = gamma, lambda = lambda, nu = nu,
#'                          zeta_mu = zeta_mu, zeta_sigma2 = zeta_sigma2)
#'
#' # Example 5
#'
#' I <- 20
#' J <- 10
#' K <- 1
#' M <- 2
#' N <- 2
#' omega_mu <- matrix(data = c(2.50, -2.00, 0.50, 0.00), nrow = 1, ncol = M * N)
#' omega_sigma2 <- diag(x = c(0.90, 0.70, 0.30, 0.10), nrow = M * N)
#' contrast_codes <- cbind(1, contr.poly(n = J))[, 1:N]
#' gamma <- matrix(data = 0, nrow = J * M, ncol = M * N)
#' for(j in 1:J) {
#'   for(m in 1:M) {
#'     gamma[(m + M * (j - 1)), (((m - 1) * N + 1):((m - 1) * N + N))] <-
#'       contrast_codes[j, ]
#'   }
#' }
#' key <- rbinom(n = I * J, size = 1, prob = .7) + 1
#' measure_weights <-
#'   matrix(data = c(0.5, -1.0, 0.5, 1.0), nrow = 2, ncol = M, byrow = TRUE)
#' lambda <- matrix(data = 0, nrow = I * J, ncol = J * M)
#' for(j in 1:J) {
#'   lambda[(1 + (j - 1) * I):(j * I), (1 + (j - 1) * M):(j * M)] <-
#'     measure_weights[key, ][(1 + (j - 1) * I):(j * I), ]
#' }
#' zeta_mu <- matrix(data = rep(x = 0, times = M * J), nrow = 1, ncol = J * M)
#' zeta_sigma2 <- diag(x = 0.2, nrow = J * M, ncol = J * M)
#' nu_mu <- matrix(data = c(0.00), nrow = 1, ncol = 1)
#' nu_sigma2 <- matrix(data = c(0.20), nrow = 1, ncol = 1)
#' set.seed(624)
#' ex5 <- dich_response_sim(I = I, J = J, K = K, M = M, N = N,
#'                          omega_mu = omega_mu, omega_sigma2 = omega_sigma2,
#'                          gamma = gamma, lambda = lambda, nu_mu = nu_mu,
#'                          nu_sigma2 = nu_sigma2, zeta_mu = zeta_mu,
#'                          zeta_sigma2 = zeta_sigma2, key = key)
#'
#' @export dich_response_sim
#-------------------------------------------------------------------------------

dich_response_sim <- function(I = NULL, J = NULL, K = NULL, M = NULL, N = NULL,
                              omega = NULL, omega_mu = NULL,
                              omega_sigma2 = NULL, gamma = NULL,
                              lambda = NULL, lambda_mu = NULL,
                              lambda_sigma2 = NULL, nu = NULL, nu_mu = NULL,
                              nu_sigma2 = NULL, zeta = NULL, zeta_mu = NULL,
                              zeta_sigma2 = NULL, kappa = NULL, key = NULL,
                              link  = "probit") {
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
  if (is.null(x = nu)) {
    nu <- matrix(
      data = MASS::mvrnorm(
        n = I * J,
        mu = nu_mu,
        Sigma = nu_sigma2
      ),
      nrow = I * J,
      ncol = 1,
      byrow = TRUE
    )
  }
  if (is.null(x = omega)) {
    omega <-  matrix(
      data = MASS::mvrnorm(
        n = K,
        mu = omega_mu,
        Sigma = omega_sigma2
      ),
      nrow = K,
      ncol = N * M,
      byrow = FALSE
    )
  }
  if (is.null(x = lambda)) {
    lambda <-  matrix(
      data = MASS::mvrnorm(
        n = I,
        mu = lambda_mu,
        Sigma = lambda_sigma2
      ),
      nrow = I * J,
      ncol = J * M,
      byrow = FALSE
    )
  }
  if (is.null(x = zeta)) {
    zeta <- matrix(
      data = MASS::mvrnorm(
        n = K,
        mu = zeta_mu,
        Sigma = zeta_sigma2 * diag(x = 1, nrow = M * J)
      ),
      nrow = K,
      ncol = J * M,
      byrow = FALSE
    )
  }
  kappa_mat <- if (is.null(x = kappa)) {
    array(data = 0, dim = c(K, I * J))
  } else {
    array(data = 1, dim = c(K, 1)) %*% t(kappa)
  }
  ystar <- array(data = 1, dim = c(K, 1)) %*% t(nu) +
    omega %*% t(gamma) %*% t(lambda) + zeta %*% t(lambda)
  p <- if (link == "logit") {
    kappa_mat + (1 - kappa_mat) * plogis(q = ystar)
  } else if (link == "probit") {
    kappa_mat + (1 - kappa_mat) * pnorm(q = ystar)
  }
  y <- apply(data.matrix(
    p > array(data = runif(n = K * I * J), dim = c(K, I * J))
  ), 1:2, as.numeric)
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
  simdat <- list(
    "y" = y,
    "ystar" = ystar,
    "omega" = omega,
    "omega_mu" = omega_mu,
    "omega_sigma2" = omega_sigma2,
    "gamma" = gamma,
    "lambda" = lambda,
    "lambda_mu" = lambda_mu,
    "lambda_sigma2" = lambda_sigma2,
    "nu" = nu,
    "nu_mu" = nu_mu,
    "nu_sigma2" = nu_sigma2,
    "zeta" = zeta,
    "zeta_mu" = zeta_mu,
    "zeta_sigma2" = zeta_sigma2,
    "kappa" = kappa,
    "condition" = condition, "key" = key
  )
  return(simdat)
}
