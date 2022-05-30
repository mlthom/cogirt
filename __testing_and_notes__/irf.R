#-------------------------------------------------------------------------------
#' Explore Conditional Item Response Functions
#-------------------------------------------------------------------------------

rm(list = ls())

# Version 1: Simulation until I get actual IRF code worked out. ----------------

xlim <- c(-5, 5)
ylim <- c(0, 1)
plot(NULL, xlim = xlim, ylim = ylim)

for(theta in seq(xlim[1], xlim[2], .1)) {

  I <- 1
  J <- 5
  K <- 500
  M <- 2
  N <- 3
nu_mu <- 0
nu_sigma2 <- 0.0
omega_mu <- matrix(data = c(theta, -0.5, 0.0, 0.0, 0.0, 0.0), nrow = 1,
                   ncol = M * N)
omega_sigma2 <- diag(x = c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00), nrow = M * N)
zeta_mu <- matrix(data = rep(x = 0, times = M * J), nrow = 1, ncol = J * M)
zeta_sigma2 <- diag(x = 0.2, nrow = J * M, ncol = J * M)
item_type <- rep(1, I * J) #1 = targ 2 = foil
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

tmp <- dich_response_sim(I = I, J = J, K = K, M = M, N = N, nu_mu = nu_mu,
                                  nu_sigma2 = nu_sigma2, lambda = lambda, gamma = gamma,
                                  omega_mu = omega_mu, omega_sigma2 = omega_sigma2,
                                  zeta_mu = zeta_mu, zeta_sigma2 = zeta_sigma2)



tmp <- dich_response_model(y = NULL, nu = tmp$nu, lambda = tmp$lambda,
                    gamma = tmp$gamma, omega = tmp$omega, zeta = tmp$zeta,
                    link = "probit")$p

points(theta, mean(tmp[, 1]), col = rainbow(5)[1])
points(theta, mean(tmp[, 2]), col = rainbow(5)[2])
points(theta, mean(tmp[, 3]), col = rainbow(5)[3])
points(theta, mean(tmp[, 4]), col = rainbow(5)[4])
points(theta, mean(tmp[, 5]), col = rainbow(5)[5])
legend(x = xlim[2], y = ylim[1], xjust = 1, yjust = 0, legend = 1:5, rainbow(5))

}
abline(h = .5)

