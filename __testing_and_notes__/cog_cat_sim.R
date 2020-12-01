#-------------------------------------------------------------------------------
#' Cognitive Testing Using Computerized Adaptive Testing Simulation
#'
#' This function simulates adaptive testing using the cat_cog.R function.
#'
#' @param x .rda file (or list) containing all objects necessary to run rmmh.
#' @param obj_fun A function that calculates predictions and log-likelihood
#' values for the selected model (character).
#' @param int_par Intentional parameters. That is, the parameters to optimize
#' precision (scalar).
#' @param min_se Minimum standard error for each intentional parameter (vector).
#'
#' @examples
#' cog_cat_sim(rda = sdirtSS, obj_fun = dich_response_model, int_par = 1,
#' min_se = -Inf, verbose_sim = T)
#'
#' @export cog_cat_sim
#'
#-------------------------------------------------------------------------------

cog_cat_sim <- function(rda = NULL, obj_fun = NULL, int_par = NULL,
                        min_se = NULL, verbose_sim = T) {
  if(verbose_sim) {
    cat(
      "CAT Simulation Start Time",
      format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
      "\n",
      sep = " "
    )
    xlim <- c(
      rda$omega[1, int_par] - 2 * rda$omega_sigma2[int_par, int_par],
      rda$omega[1, int_par] + 2 * rda$omega_sigma2[int_par, int_par]
      )
    xseg <- (xlim[2] - xlim[1]) / 5
    ylim <- c(0, rda$omega_sigma2[int_par, int_par])
    yseg <- rda$omega_sigma2[int_par, int_par] / 5
    plot(NULL, axes = F, xlim = xlim, ylim = ylim, xlab = "Parameter", ylab = "Standard Error", main = "")
    axis(side = 1, tick = T, at = seq(xlim[1],xlim[2],xseg))
    axis(side = 2, tick = T, at = seq(ylim[1],ylim[2],yseg))
    abline(v = rda$omega[1, int_par], lwd = 2, col = "green")
    abline(v = rda$omega_mu[1, int_par], lwd = 2, col = "red")
    tmp <- rmmh(
      chains = 3,
      y = rda$y,
      obj_fun = obj_fun,
      est_omega = T,
      est_nu = T,
      est_zeta = T,
      lambda = rda$lambda,
      gamma = rda$gamma,
      omega0 = array(data = 0, dim = dim(rda$omega)),
      nu0 = array(
        data = 0,
        dim = c(ncol(rda$nu), 1)
      ),
      zeta0 = array(data = 0, dim = dim(rda$zeta)),
      omega_mu = rda$omega_mu,
      omega_sigma2 = rda$omega_sigma2,
      nu_mu = matrix(rda$nu_mu),
      nu_sigma2 = matrix(rda$nu_sigma2),
      zeta_mu = rda$zeta_mu,
      zeta_sigma2 = rda$zeta_sigma2,
      burn = 0,
      thin = 1,
      min_tune = 0,
      tune_int = Inf,
      max_tune = 0,
      niter = 1,
      verbose_rmmh = F,
      max_iter_rmmh = 200
    )
    abline(v = tmp$omega1[1, int_par], lwd = 2, col = "yellow")
  }
  rda$list <- c(sapply(X = 1:(length(rda$y) / 5), FUN = rep, 5))
  rda_sim <- rda
  rda_sim$y[which(!rda_sim$list %in% c(1))] <- NA
  se <- matrix(data = Inf, nrow = length(int_par), ncol = length(int_par))
  iter <- 0
  while (any(diag(x = se) > min_se) & any(is.na(x = rda_sim$y))) {
    iter <- iter + 1
    tmp <- cog_cat(rda = rda_sim, obj_fun = obj_fun, int_par = int_par)
    se <- matrix(
      data = solve(tmp$info1)[int_par, int_par],
      nrow = length(int_par),
      ncol = length(int_par)
    )
    rda_sim[["y"]][which(rda_sim$list == tmp$next_list)] <-
      rda[["y"]][which(rda$list == tmp$next_list)]
    if(verbose_sim) {
      cat(
        "... at iteration ",
        format(x = iter, nsmall = 0),
        " MAP omega is ",
        format(x = round(x = tmp$omega1[1, int_par], digits = 3), nsmall = 3),
        " true value is ",
        format(x = round(x = rda$omega[1, int_par], digits = 3), nsmall = 3),
        " PSD omega is ",
        format(x = round(x = se, digits = 3), nsmall = 3),
        " next list is ",
        tmp$next_list,
        " next condition is ",
        unique(x = rda$condition[which(rda$list == tmp$next_list)]),
        "\n",
        sep = ""
      )
      points(x = tmp$omega1[1, int_par], y = se, pch = 15, col = "blue", xpd = T)
      segments(
        x0 = tmp$omega1[1, int_par] - se,
        y0 = se,
        x1 = tmp$omega1[1, int_par] + se,
        y1 = se,
        col = "blue",
        xpd = T
        )
    }
  }
  return(list(
    "omega1" = tmp$omega1,
    "se" = se
  ))
}
