#-------------------------------------------------------------------------------
#' Compare RMMH and MCMH Estimators In an Adaptive Testing Scenario
#'
#' This script is used only for developmental purposes.
#'
#-------------------------------------------------------------------------------

obj_fun = dich_response_model
int_par = 1
min_se = -Inf
verbose_sim = T
max_list = 5

num_sim <- 1000
sim_res <- NULL
for(ns in 1:num_sim){

  # STEP 1: Simulate data for single subject -------------------------------------

  I = 200
  J = 5
  K = 1
  M = 2
  N = 3
  nu_mu = 0
  nu_sigma2 = 0.2
  omega_mu <- matrix(data = c(3.00, -2.00, 0.10, 0.00, 0.50, -0.10), nrow = 1,
                     ncol = M * N)
  omega_sigma2 <- diag(x = c(0.50, 0.10, 0.01, 0.25, 0.10, 0.01), nrow = M * N)
  zeta_mu <- matrix(data = rep(x = 0, times = M * J), nrow = 1, ncol = J * M)
  zeta_sigma2 <- diag(x = 0.2, nrow = J * M, ncol = J * M)
  item_type <- rbinom(n = I * J, size = 1, prob = .7) + 1
  # Equation 12 Thomas et al. (2018)
  measure_weights <-
    matrix(data = c(0.5, -1.0, 0.5, 1.0), nrow = 2, ncol = M, byrow = T)
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

  # STEP 2: Set up plot and estimate predicted score given compelte data ---------

  xlim <- c(
    rda$omega_mu[1, int_par] - 3 * rda$omega_sigma2[int_par, int_par],
    rda$omega_mu[1, int_par] + 3 * rda$omega_sigma2[int_par, int_par]
  )
  xseg <- (xlim[2] - xlim[1]) / 5
  ylim <- c(0, rda$omega_sigma2[int_par, int_par])
  yseg <- rda$omega_sigma2[int_par, int_par] / 5
  par(mfrow = c(1,1))
  plot(NULL, axes = F, xlim = xlim, ylim = ylim, xlab = "Parameter",
       ylab = "Standard Error", main = "")
  axis(side = 1, tick = T, at = seq(xlim[1],xlim[2],xseg))
  axis(side = 2, tick = T, at = seq(ylim[1],ylim[2],yseg))
  abline(v = rda$omega[1, int_par], lwd = 2, col = "green")
  abline(v = rda$omega_mu[1, int_par], lwd = 2, col = "red")
  tmp <- rmmh(chains = 3, y = rda$y, obj_fun = obj_fun, est_omega = T,
              est_nu = T, est_zeta = T, lambda = rda$lambda, gamma = rda$gamma,
              omega0 = array(data = 0, dim = dim(rda$omega)),
              nu0 = array(data = 0, dim = c(ncol(rda$nu), 1)),
              zeta0 = array(data = 0, dim = dim(rda$zeta)), omega_mu = rda$omega_mu,
              omega_sigma2 = rda$omega_sigma2, nu_mu = matrix(rda$nu_mu),
              nu_sigma2 = matrix(rda$nu_sigma2), zeta_mu = rda$zeta_mu,
              zeta_sigma2 = rda$zeta_sigma2, burn = 0, thin = 1, min_tune = 0,
              tune_int = Inf, max_tune = 0, niter = 1, verbose_rmmh = F,
              max_iter_rmmh = 200
  )
  abline(v = tmp$omega1[1, int_par], lwd = 2, col = "purple")

  text(x = rda$omega[1, int_par], y = ylim[2], labels = "True Value", xpd = T,
       srt = 45)
  text(x = rda$omega_mu[1, int_par], y = ylim[2], labels = "Prior Value", xpd = T,
       srt = 45)
  text(x = tmp$omega1[1, int_par], y = ylim[2], labels = "Complee Data Value",
       xpd = T, srt = 45)

  # Step 3: Estimate trait using adaptive testing --------------------------------

  rda$list <- c(sapply(X = 1:(length(rda$y) / 20), FUN = rep, 20))
  rda_sim <- rda
  rda_sim$y[which(!rda_sim$list %in% c(2))] <- NA
  adap_se <- matrix(data = Inf, nrow = length(int_par), ncol = length(int_par))
  iter <- 0
  while (any(diag(x = adap_se) > min_se) & any(is.na(x = rda_sim$y)) & iter < max_list) {
    iter <- iter + 1
    adap_res <- cog_cat(rda = rda_sim, obj_fun = obj_fun, int_par = int_par)
    adap_se <- matrix(
      data = solve(adap_res$info1)[int_par, int_par],
      nrow = length(int_par),
      ncol = length(int_par)
    )
    rda_sim[["y"]][which(rda_sim$list == adap_res$next_list)] <-
      rda[["y"]][which(rda$list == adap_res$next_list)]
    if(verbose_sim) {
      cat(
        "... at iteration ",
        format(x = iter, nsmall = 0),
        " MAP omega is ",
        format(x = round(x = adap_res$omega1[1, int_par], digits = 3), nsmall = 3),
        " true value is ",
        format(x = round(x = rda$omega[1, int_par], digits = 3), nsmall = 3),
        " PSD omega is ",
        format(x = round(x = adap_se, digits = 3), nsmall = 3),
        " next list is ",
        adap_res$next_list,
        " next condition is ",
        unique(x = rda$condition[which(rda$list == adap_res$next_list)]),
        " percet correct is ",
        format(x = round(x = mean(rda[["y"]][which(rda$list == adap_res$next_list)]), digits = 2), nsmall = 2),
        "\n",
        sep = ""
      )
      points(x = adap_res$omega1[1, int_par], y = adap_se, pch = 15, col = "blue", xpd = T)
      segments(
        x0 = adap_res$omega1[1, int_par] - adap_se,
        y0 = adap_se,
        x1 = adap_res$omega1[1, int_par] + adap_se,
        y1 = adap_se,
        col = "blue",
        xpd = T
      )
    }
  }

  # Step 4: Estimate trait using non adaptive testing ----------------------------

  rda$list <- c(sapply(X = 1:(length(rda$y) / 5), FUN = rep, 5))
  rda_sim <- rda
  rda_sim$y[which(!rda_sim$list %in% c(2))] <- NA
  non_adap_se <- matrix(data = Inf, nrow = length(int_par), ncol = length(int_par))
  iter <- 0
  next_condition <- 1
  while (any(diag(x = non_adap_se) > min_se) & any(is.na(x = rda_sim$y)) & iter < max_list) {
    iter <- iter + 1
    if(next_condition < 5){
      next_condition <- next_condition + 1
    } else {
      next_condition <- 1
    }
    next_list <- rda_sim$list[min(which(rda_sim$condition == next_condition & is.na(rda_sim$y)))]
    non_adap_res <- cog_cat(rda = rda_sim, obj_fun = obj_fun, int_par = int_par)
    non_adap_se <- matrix(
      data = solve(non_adap_res$info1)[int_par, int_par],
      nrow = length(int_par),
      ncol = length(int_par)
    )
    rda_sim[["y"]][which(rda_sim$list == next_list)] <-
      rda[["y"]][which(rda$list == next_list)]
    if(verbose_sim) {
      cat(
        "... at iteration ",
        format(x = iter, nsmall = 0),
        " MAP omega is ",
        format(x = round(x = non_adap_res$omega1[1, int_par], digits = 3), nsmall = 3),
        " true value is ",
        format(x = round(x = rda$omega[1, int_par], digits = 3), nsmall = 3),
        " PSD omega is ",
        format(x = round(x = non_adap_se, digits = 3), nsmall = 3),
        " next list is ",
        next_list,
        " next condition is ",
        next_condition,
        " percet correct is ",
        format(x = round(x = mean(rda[["y"]][which(rda$list == next_list)]), digits = 2), nsmall = 2),
        "\n",
        sep = ""
      )
      points(x = non_adap_res$omega1[1, int_par], y = non_adap_se, pch = 15, col = "red", xpd = T)
      segments(
        x0 = non_adap_res$omega1[1, int_par] - non_adap_se,
        y0 = non_adap_se,
        x1 = non_adap_res$omega1[1, int_par] + non_adap_se,
        y1 = non_adap_se,
        col = "red",
        xpd = T
      )
    }
  }

  sim_res <- rbind(sim_res, c(tmp$omega1[1, int_par], rda$omega[1], adap_res$omega[1], non_adap_res$omega[1], adap_se, non_adap_se))
  print(ns)

}

sim_res <- data.frame(sim_res)
colnames(sim_res) <- c("comp_data_omega", "true_omega", "adapt_omega", "non_adapt_omega", "adapt_se", "non_adapt_se")

par(mfrow = c(1,1))
plot(sim_res$true_omega, sim_res$adapt_omega, col = "blue", pch = 1)
points(sim_res$true_omega, sim_res$non_adapt_omega, col = "red", pch = 2)
cor(sim_res$true_omega, sim_res$adapt_omega)
cor(sim_res$true_omega, sim_res$non_adapt_omega)
abline(lm(sim_res$adapt_omega ~ sim_res$true_omega), col = "blue")
abline(lm(sim_res$non_adapt_omega ~ sim_res$true_omega), col = "red")

par(mfrow = c(1,1))
plot(sim_res$true_omega, sim_res$adapt_omega - sim_res$true_omega, col = "blue", pch = 1)
points(sim_res$true_omega, sim_res$non_adapt_omega - sim_res$true_omega, col = "red", pch = 2)
abline(v = rda$omega_mu[1])
cor(sim_res$true_omega, sim_res$adapt_omega - sim_res$true_omega)
cor(sim_res$true_omega, sim_res$non_adapt_omega - sim_res$true_omega)
abline(lm(I(sim_res$adapt_omega) - sim_res$true_omega ~ sim_res$true_omega), col = "blue")
abline(lm(I(sim_res$non_adapt_omega) - sim_res$true_omega ~ sim_res$true_omega), col = "red")

par(mfrow = c(1,2))
hist(sim_res$adapt_se)
hist(sim_res$non_adapt_se)
mean(sim_res$adapt_se)
mean(sim_res$non_adapt_se)

par(mfrow = c(1,1))
plot(sim_res$true_omega, sim_res$comp_data_omega, col = "blue", pch = 1)
cor(sim_res$true_omega, sim_res$comp_data_omega)

