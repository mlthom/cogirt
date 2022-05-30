#-------------------------------------------------------------------------------
#' Compare RMMH and MCMH Estimators
#'
#' This script is used only for developmental purposes. It simulates data for a
#' single subject, and then estimates the omega parameters using the CatCOG
#' mcmh and rmmh functions. In then plots the results.
#'
#-------------------------------------------------------------------------------

# STEP 1: Simulate Data --------------------------------------------------------

if (1 == 0){

  # 6 parameters

  int_par_1 <- 1
  int_par_2 <- 4

  I = 4
  J = 5
  K = 1
  M = 2
  N = 3
  nu_mu = 0
  nu_sigma2 = 0.2
  omega_mu <- matrix(data = c(2.30, -1.90, 0.50, 0.50, 0.40, -0.10), nrow = 1,
                     ncol = M * N)
  omega_sigma2 <- diag(x = c(0.50, 0.10, 0.01, 0.25, 0.10, 0.01), nrow = M * N)
  zeta_mu <- matrix(data = rep(x = 0, times = M * J), nrow = 1, ncol = J * M)
  zeta_sigma2 <- diag(x = 0.2, nrow = J * M, ncol = J * M)
  item_type <- rbinom(n = I * J, size = 1, prob = .7) + 1
  # Equation 12 Thomas et al. (2018)
  measure_weights <-
    matrix(data = c(0.5, -1.0, 0.5, 1.0), nrow = 2, ncol = M, byrow = TRUE)
  lambda <- matrix(data = 0, nrow = I * J, ncol = J * M)
  for(j in 1:J) {
    lambda[(1 + (j - 1) * I):(j * I), (1 + (j - 1) * M):(j * M)] <-
      measure_weights[item_type, ][(1 + (j - 1) * I):(j * I), ]
  }
  contrast_codes <- cbind(1, contr.poly(n = J))[, 1:N]
  gamma <- matrix(data = 0, nrow = J * M, ncol = M * N)
  for(j in 1:J) {
    for(m in 1:M) {
      gamma[(m + M * (j - 1)), (((m - 1) * N + 1):((m - 1) * N + N))] <-
        contrast_codes[j, ]
    }
  }

  rda <- dich_response_sim(I = I, J = J, K = K, M = M, N = N,
                           nu_mu = nu_mu, nu_sigma2 = nu_sigma2, lambda = lambda,
                           gamma = gamma, omega_mu = omega_mu,
                           omega_sigma2 = omega_sigma2, zeta_mu = zeta_mu,
                           zeta_sigma2 = zeta_sigma2)
} else if (1 == 0) {

  # 4 parameters

  int_par_1 <- 1
  int_par_2 <- 3

  I = 20
  J = 5
  K = 1
  M = 2
  N = 2
  nu_mu = 0
  nu_sigma2 = 0.2
  omega_mu <- matrix(data = c(3.00, -2.00, 0.00, 0.50), nrow = 1,
                     ncol = M * N)
  omega_sigma2 <- diag(x = c(0.50, 0.10, 0.25, 0.10), nrow = M * N)
  zeta_mu <- matrix(data = rep(x = 0, times = M * J), nrow = 1, ncol = J * M)
  zeta_sigma2 <- diag(x = 0.2, nrow = J * M, ncol = J * M)
  item_type <- rbinom(n = I * J, size = 1, prob = .7) + 1
  # Equation 12 Thomas et al. (2018)
  measure_weights <-
    matrix(data = c(0.5, -1.0, 0.5, 1.0), nrow = 2, ncol = M, byrow = TRUE)
  lambda <- matrix(data = 0, nrow = I * J, ncol = J * M)
  for(j in 1:J) {
    lambda[(1 + (j - 1) * I):(j * I), (1 + (j - 1) * M):(j * M)] <-
      measure_weights[item_type, ][(1 + (j - 1) * I):(j * I), ]
  }
  contrast_codes <- cbind(1, contr.poly(n = J))[, 1:N]
  gamma <- matrix(data = 0, nrow = J * M, ncol = M * N)
  for(j in 1:J) {
    for(m in 1:M) {
      gamma[(m + M * (j - 1)), (((m - 1) * N + 1):((m - 1) * N + N))] <-
        contrast_codes[j, ]
    }
  }

  rda <- dich_response_sim(I = I, J = J, K = K, M = M, N = N,
                           nu_mu = nu_mu, nu_sigma2 = nu_sigma2, lambda = lambda,
                           gamma = gamma, omega_mu = omega_mu,
                           omega_sigma2 = omega_sigma2, zeta_mu = zeta_mu,
                           zeta_sigma2 = zeta_sigma2)
} else if (1 == 0) {

  # 2 parameters (note that mean and variance are changed here)

  int_par_1 <- 1
  int_par_2 <- 2

  I = 20
  J = 5
  K = 1
  M = 2
  N = 2
  nu_mu = 0
  nu_sigma2 = 0.2
  omega_mu <- matrix(data = c(0.00, 0.00), nrow = 1,
                     ncol = M * N)
  omega_sigma2 <- diag(x = c(1.00, 1.00), nrow = M * N)
  zeta_mu <- matrix(data = rep(x = 0, times = M * J), nrow = 1, ncol = J * M)
  zeta_sigma2 <- diag(x = 0.2, nrow = J * M, ncol = J * M)
  item_type <- rbinom(n = I * J, size = 1, prob = .7) + 1
  # Equation 12 Thomas et al. (2018)
  measure_weights <-
    matrix(data = c(0.5, -1.0, 0.5, 1.0), nrow = 2, ncol = M, byrow = TRUE)
  lambda <- matrix(data = 0, nrow = I * J, ncol = J * M)
  for(j in 1:J) {
    lambda[(1 + (j - 1) * I):(j * I), (1 + (j - 1) * M):(j * M)] <-
      measure_weights[item_type, ][(1 + (j - 1) * I):(j * I), ]
  }
  contrast_codes <- cbind(1, contr.poly(n = J))[, 1:N]
  gamma <- matrix(data = 0, nrow = J * M, ncol = M * N)
  for(j in 1:J) {
    for(m in 1:M) {
      gamma[(m + M * (j - 1)), (((m - 1) * N + 1):((m - 1) * N + N))] <-
        contrast_codes[j, ]
    }
  }

  rda <- dich_response_sim(I = I, J = J, K = K, M = M, N = N,
                           nu_mu = nu_mu, nu_sigma2 = nu_sigma2, lambda = lambda,
                           gamma = gamma, omega_mu = omega_mu,
                           omega_sigma2 = omega_sigma2, zeta_mu = zeta_mu,
                           zeta_sigma2 = zeta_sigma2)
} else if( 1 == 1 ) {

  # Guessing model

  int_par_1 <- 1
  int_par_2 <- 2

  I <- 5
  J <- 5
  K <- 1
  M <- 1
  N <- 2
  nu_mu <- 0
  nu_sigma2 <- 0.2
  omega_mu <- matrix(data = c(-1.00, 1.00), nrow = 1, ncol = M * N)
  omega_sigma2 <- diag(x = c(1.00, 0.01), nrow = M * N)
  zeta_mu <- matrix(data = rep(x = 0, times = M * J), nrow = 1, ncol = J * M)
  zeta_sigma2 <- diag(x = 0.2, nrow = J * M, ncol = J * M)
  lambda <- matrix(data = 0, nrow = I * J, ncol = J * M)
  for (j in 1:J) {
    lambda[(1 + (j - 1) * I):(j * I), j] <- 1
  }
  contrast_codes <- cbind(1, contr.poly(n = J))[, 1:N]
  gamma <- matrix(data = 0, nrow = J * M, ncol = M * N)
  for(j in 1:J) {
    for(m in 1:M) {
      gamma[(m + M * (j - 1)), (((m - 1) * N + 1):((m - 1) * N + N))] <-
        contrast_codes[j, ]
    }
  }
  kappa <- array(data = 0.5, dim = c(K, I * J))

  rda <- dich_response_sim(I = I, J = J, K = K, M = M, N = N,
                                   nu_mu = nu_mu, nu_sigma2 = nu_sigma2,
                                   lambda = lambda, kappa = kappa,
                                   gamma = gamma, omega_mu = omega_mu,
                                   omega_sigma2 = omega_sigma2,
                                   zeta_mu = zeta_mu,
                                   zeta_sigma2 = zeta_sigma2, item_type = NULL)
}

# STEP 2: MCMH Estimation ------------------------------------------------------

mcmh_time <- system.time(mcmh_est <- mcmh_mc(
  chains = 3,
  y = rda$y,
  obj_fun = dich_response_model,
  est_omega = TRUE,
  est_nu = TRUE,
  est_zeta = TRUE,
  lambda = rda$lambda,
  kappa = rda$kappa,
  gamma = rda$gamma,
  omega0 = rda$omega_mu,
  nu0 = array(data = 0, dim = c(ncol(rda$nu), 1)),
  zeta0 = array(data = 0, dim = dim(rda$zeta)),
  omega_mu = rda$omega_mu,
  omega_sigma2 = rda$omega_sigma2,
  nu_mu = matrix(rda$nu_mu),
  nu_sigma2 = matrix(rda$nu_sigma2),
  zeta_mu = rda$zeta_mu,
  zeta_sigma2 = rda$zeta_sigma2,
  burn = 0,
  thin = 10,
  min_tune = 50,
  tune_int = 50,
  max_tune = 500,
  niter = 1000
))

# STEP 3: RMMH Estimation ------------------------------------------------------

rmmh_time <- system.time(rmmh_est <- rmmh(
  chains = 1,
  y = rda$y,
  obj_fun = dich_response_model,
  est_omega = TRUE,
  est_nu = TRUE,
  est_zeta = TRUE,
  lambda = rda$lambda,
  kappa = rda$kappa,
  gamma = rda$gamma,
  omega0 = rda$omega_mu,
  nu0 = array(data = 0, dim = c(ncol(rda$nu), 1)),
  zeta0 = array(data = 0, dim = dim(rda$zeta)),
  omega_mu = rda$omega_mu, omega_sigma2 = rda$omega_sigma2,
  nu_mu = matrix(rda$nu_mu), nu_sigma2 = matrix(rda$nu_sigma2),
  zeta_mu = rda$zeta_mu,
  zeta_sigma2 = rda$zeta_sigma2,
  burn=5,
  thin=1,
  min_tune=0,
  tune_int=0,
  max_tune = 0,
  niter = 6,
  verbose_rmmh = TRUE,
  max_iter_rmmh = 200
))

# STEP 4: NR Estimation --------------------------------------------------------

nr_time <- system.time(nr_est <- nr(rda = rda, verbose = F))

# STEP 5: Plot Estimates -------------------------------------------------------

xlim <- c(
  rda$omega_mu[1, int_par_1] - 4 * rda$omega_sigma2[int_par_1, int_par_1],
  rda$omega_mu[1, int_par_1] + 4 * rda$omega_sigma2[int_par_1, int_par_1]
)
xseg <- (xlim[2] - xlim[1]) / 5

ylim <- c(
  rda$omega_mu[1, int_par_2] - 10 * rda$omega_sigma2[int_par_2, int_par_2],
  rda$omega_mu[1, int_par_2] + 10 * rda$omega_sigma2[int_par_2, int_par_2]
)
yseg <- (ylim[2] - ylim[1]) / 5

plot(NULL, axes = F, xlim = xlim, ylim = ylim, xlab = "Intentional Parameter 1",
     ylab = "Intentional Parameter 2", main = "")
rect(xleft = -1000, xright = 1000, ybottom = -1000, ytop = 1000, border = NA,
     col = "grey90")
abline(v = seq(xlim[1], xlim[2], xseg), lwd = 2, col = "white")
abline(h = seq(ylim[1], ylim[2], yseg), lwd = 2, col = "white")
axis(side = 1, tick = F, at = seq(xlim[1], xlim[2], xseg))
axis(side = 2, tick = F, at = seq(ylim[1], ylim[2], yseg))

points(x = rda$omega[int_par_1], y = rda$omega[int_par_2], col = rgb(1,0,0,.5),
       pch = 20, cex = 4, xpd = TRUE)
points(x = mcmh_est$omegaEAP[int_par_1], y = mcmh_est$omegaEAP[int_par_2],
       col = rgb(0,1,0,.5), pch = 20, cex = 4, xpd = TRUE)
points(x = rmmh_est$omega1[int_par_1], y = rmmh_est$omega1[int_par_2],
       col = rgb(0,0,1,.5), pch = 20, cex = 4, xpd = TRUE)
points(x = nr_est[int_par_1], y = nr_est[int_par_2],
       col = rgb(1,0,1,.5), pch = 20, cex = 4, xpd = TRUE)

legend(x = xlim[1], y = ylim[2], yjust = -0.1, xpd = TRUE, y.intersp = .8, legend = c(
  "True",
  paste("MC-MH ", round(mcmh_time[3], 2), "\"", sep = ""),
  paste("RM-MH ", round(rmmh_time[3], 2), "\"", sep = ""),
  paste("Newton Raphson ", round(nr_time[3], 2), "\"", sep = "")
),
fill = c(rgb(1,0,0,.5), rgb(0,1,0,.5), rgb(0,0,1,.5), rgb(1,0,1,.5)), bty = "n",
cex = 1.00, border = c("red", "green", "blue", "purple"))





