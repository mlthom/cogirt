#-------------------------------------------------------------------------------
#' MCMH Parameter Estimates for Single Chain
#'
#' This function calculates MCMH parameter estimates for a single chain. The
#' MCMH method implemented follows Patz and Junker (1999). The approach to
#' tuning the scale and covariance matrix follows BDA and SAS 9.2 User Guide,
#' 2nd Ed. "The MCMC Procedure: Tuning the Proposal Distribution".
#'
#' @param y Matrix of item responses (K by IJ).
#' @param obj_fun A function that calculates predictions and log-likelihood
#' values for the selected model (character).
#' @param est_omega Determines whether omega is estimated (logical).
#' @param est_nu Determines whether nu is estimated (logical).
#' @param est_zeta Determines whether zeta is estimated (logical).
#' @param lambda Matrix of item structure parameters (IJ by JM).
#' @param gamma Matrix of experimental structure parameters (JM by MN).
#' @param omega0 Starting values for omega.
#' @param nu0 Starting values for nu.
#' @param zeta0 Starting values for zeta.
#' @param omega_mu Vector of means prior for omega (1 by MN).
#' @param omega_sigma Covariance matrix prior for omega (MN by MN).
#' @param zeta_mu Vector of means prior for zeta (1 by JM).
#' @param zeta_sigma@ Covariance matrix prior for zeta (JM by JM).
#' @param nu_mu Prior mean for nu (scalar).
#' @param nu_sigma@ Prior variance for nu (scalar).
#' @param burn Number of iterations at the beginning of an MCMC run to discard
#' (scalar).
#' @param thin Determines every nth observation retained (scalar).
#' @param min_tune Determines when tunning begins (scalar).
#' @param tune_int MCMH tuning interval (scalar).
#' @param max_tune Determines when tunning ends (scalar).
#' @param niter Number of iterations of the MCMH sampler.
#' @param weight Determines the weight of old versus new covariance matrix.
#' @param verbose_mcmh Print progress of MCMH sampler.
#'
#' @return List with elements omega_draws (list of (niter - burn) / thin draws
#' for K by MN omega matrix), nu_draws (list of (niter - burn) / thin draws
#' for K by IJ nu matrix), zeta_draws (list of (niter - burn) / thin draws for K
#' by JM zeta matrix), cand_o_var (list of K final MN by MN candidate proposal
#' covariance matrices for omega for each examinee), cand_n_var (list of IJ
#' final scalar candidate proposal variances for nu for all items), cand_z_var
#' (list of final JM by JM candidate proposal covariance matrices for zeta for
#' all examinees)
#'
#' @references
#'
#' Patz, R. J., & Junker, B. W. (1999). A Straightforward Approach to Markov
#' Chain Monte Carlo Methods for Item Response Models. Journal of Educational
#' and Behavioral Statistics, 24(2), 146.
#'
#' @examples
#'mcmh_sc(y = sdirt$y, obj_fun = dich_response_model, est_omega = T,
#'     est_nu = T, est_zeta = T, lambda = sdirt$lambda, gamma = sdirt$gamma,
#'     omega0 = array(data = 0, dim = dim(sdirt$omega)),
#'     nu0 = array(data = 0, dim = c(ncol(sdirt$nu), 1)),
#'     zeta0 = array(data = 0, dim = dim(sdirt$zeta)),
#'     omega_mu = sdirt$omega_mu, omega_sigma2 = sdirt$omega_sigma2,
#'     nu_mu = matrix(sdirt$nu_mu), nu_sigma2 = matrix(sdirt$nu_sigma2),
#'     zeta_mu = sdirt$zeta_mu, zeta_sigma2 = sdirt$zeta_sigma2,
#'     burn = 0, thin = 10, min_tune = 50, tune_int = 50, max_tune = 1000,
#'     niter = 2000, weight = 1/1, verbose_mcmh = T)
#'
#' @export mcmh_sc
#-------------------------------------------------------------------------------

mcmh_sc <- function(
  y = y, obj_fun = NULL, est_omega = T, est_nu = T, est_zeta = T, lambda = NULL,
  gamma = NULL, omega0 = NULL, nu0 = NULL, zeta0 = NULL, omega_mu = NULL,
  omega_sigma2 = NULL, nu_mu = NULL, nu_sigma2 = NULL, zeta_mu = NULL,
  zeta_sigma2 = NULL, burn = NULL, thin = NULL, min_tune = NULL,
  tune_int = NULL, max_tune = NULL, niter = NULL, weight = 1,
  verbose_mcmh = F
  ) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package \"MASS\" needed for the mcmh_sc function to work. Please
         install.",
         call. = FALSE)
  }
  if (!requireNamespace("abind", quietly = TRUE)) {
    stop("Package \"abind\" needed for the mcmh_sc function to work. Please
         install.",
         call. = FALSE)
  }
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package \"mvtnorm\" needed for the mcmh_sc function to work. Please
         install.",
         call. = FALSE)
  }

  if (verbose_mcmh) {
    cat(
      "\n",
      "MCMC LOG",
      "Start Time",
      format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
      "\n\n",
      sep = " "
    )
  }
  #Initialize objects to save recent draws for tuning
  recent_omega_draws <- array(
    data = t(omega0),
    dim = c(1, ncol(omega0), nrow(omega0))
  )
  recent_nu_draws <- array(
    data = t(omega0),
    dim = c(1, 1, nrow(nu0))
  )
  recent_zeta_draws <- array(
    data = t(zeta0),
    dim = c(1, ncol(zeta0), nrow(zeta0))
  )
  # Initialize objects to save final draws
  omega_draws <- NULL
  nu_draws <- NULL
  zeta_draws <- NULL
  # Candidate scale
  o_scale <- 2.38 / sqrt(x = ncol(x = omega0))
  n.scale <- 2.38 / sqrt(x = ncol(x = nu0))
  z.scale <- 2.38 / sqrt(x = ncol(x = zeta0))
  # Candidate proposal variance (see BDA3 p. 296)
  o.prop <- array(
    data = diag(.75, ncol(omega0)),
    dim = c(ncol(omega0), ncol(omega0), nrow(omega0))
  )
  cand_o_var <- o.prop * o_scale
  n.prop <- array(
    data = diag(.25, 1),
    dim = c(1, 1, nrow(nu0))
  )
  cand.n.var <- n.prop * n.scale
  z.prop <- array(
    data = diag(.05, ncol(zeta0)),
    dim = c(ncol(zeta0), ncol(zeta0), nrow(zeta0))
  )
  cand_z_var <- z.prop * z.scale
  for (i in 1:niter) {
    iter_num <- i
    # Draw omegas
    if (est_omega) {
      omega1 <- omega0 + matrix(
        data = sapply(
          X = seq_len(length.out = nrow(omega0)),
          FUN = function(i) {
            MASS::mvrnorm(
              n = 1,
              mu = rep(0, ncol(omega0)),
              Sigma = cand_o_var[, , i])
          }
        ),
        nrow = nrow(omega0),
        ncol = ncol(omega0),
        byrow = T
      )
      p0 <- obj_fun(y = y, nu = array(data = nu0, dim = dim(y)),
                    lambda = lambda, gamma = gamma, omega = omega0,
                    zeta = zeta0)$p
      tmp0 <- rowSums(x = log((p0^y) * (1 - p0) ^ (1 - y)), na.rm = T)
      ll0 <- tmp0 + log(x = mvtnorm::dmvnorm(
        x = omega0,
        mean = omega_mu,
        sigma = omega_sigma2
      ))
      p1 <- obj_fun(y = y, nu = array(data = nu0, dim = dim(y)),
                    lambda = lambda, gamma = gamma, omega = omega1,
                    zeta = zeta0)$p
      tmp1 <- rowSums(x = log((p1^y) * (1 - p1) ^ (1 - y)), na.rm = T)
      ll1 <- tmp1 + log(x = mvtnorm::dmvnorm(
        x = omega1,
        mean = omega_mu,
        sigma = omega_sigma2
      ))
      # Accept draws
      accept <- ll1 - ll0
      accept <- ifelse(test = accept > 0, yes = 0, no = accept)
      accept <- ifelse(
        test = runif(nrow(omega0)) < exp(x = accept),
        yes = 1,
        no = 0
      )
      omega1[which(x = accept != 1), ] <- omega0[which(x = accept != 1), ]
      if (iter_num < tune_int) {
        recent_omega_draws <- abind::abind(
          recent_omega_draws,
          array(data = t(omega1), dim = c(1, ncol(omega0), nrow(omega0))),
          along = 1
        )
      } else {
        recent_omega_draws <- abind::abind(
          array(
            recent_omega_draws[-1, , ],
            dim(x = recent_omega_draws) - c(1, 0, 0)
            ),
          array(data = t(x = omega1), dim = c(1, ncol(omega0), nrow(omega0))),
          along = 1
        )
      }
      # Adjust tuning rate and candidate covariance matrix
      omegarate <- apply(
        X = recent_omega_draws,
        MARGIN = c(3),
        FUN = function(x) {
          mean(x = !duplicated(x = x))
        }
      )
      if (iter_num %in% seq(from = min_tune, to = max_tune, by = tune_int)) {
        o_scale <- (o_scale * qnorm(p = .23 / 2)) /
          qnorm(p = (omegarate * .997 + .001) / 2)
        cand_o_var <- (
          (
            (iter_num / (max_tune + weight * max_tune)) *
              array(
                data = apply(
                  X = array(
                    data = recent_omega_draws[
                      seq(
                        from = 1,
                        to = tune_int,
                        by = ceiling(x = tune_int / 10)
                      ), ,
                      ],
                    dim = c(
                      length(x = seq(
                        from = 1,
                        to = tune_int,
                        by = ceiling(x = tune_int / 10))
                      ),
                      ncol(omega0),
                      nrow(omega0)
                    )
                  ),
                  MARGIN = c(3),
                  FUN = var
                ),
                dim = c(ncol(omega0), ncol(omega0), nrow(omega0))
              )
          ) +
            ((1 - (iter_num / (max_tune + weight * max_tune)))) *
            array(
              data = unlist(lapply(X = o_scale, function(x) {
                diag(x = x, ncol(omega0))
              }
              )
              ),
              dim = c(ncol(omega0), ncol(omega0), nrow(omega0))
            )
        )
      }
      if (iter_num > burn & iter_num %in% seq(
        from = burn + 1,
        to = niter,
        by = thin
      )
      ) {
        omega_draws <- abind::abind(
          omega_draws,
          array(data = t(omega1), dim = c(1, ncol(omega0), nrow(omega0))),
          along = 1)
      }
    } else {
      omega1 <- omega0
    }
    # Draw nus
    if (est_nu) {
      nu1 <- nu0 + matrix(
        data = sapply(
          X = seq_len(length.out = nrow(nu0)),
          FUN = function(j) {
            MASS::mvrnorm(
              n = 1,
              mu = rep(0, 1),
              Sigma = cand.n.var[, , j]
            )
          }
        ), nrow = nrow(nu0),
        ncol = ncol(nu0),
        byrow = T
      )
      p0 <- obj_fun(y = y, nu = array(data = nu0, dim = dim(y)),
                    lambda = lambda, gamma = gamma, omega = omega0,
                    zeta = zeta0)$p
      tmp0 <- colSums(x = log((p0^y) * (1 - p0) ^ (1 - y)), na.rm = T)
      ll0 <- tmp0 + log(x = mvtnorm::dmvnorm(
        x = nu0,
        mean = nu_mu,
        sigma = nu_sigma2
      ))
      p1 <- obj_fun(y = y, nu = array(data = nu1, dim = dim(y)),
                    lambda = lambda, gamma = gamma, omega = omega0,
                    zeta = zeta0)$p
      tmp1 <- colSums(x = log((p1 ^ y) * (1 - p1) ^ (1 - y)), na.rm = T)
      ll1 <- tmp1 + log(x = mvtnorm::dmvnorm(
        x = nu1,
        mean = nu_mu,
        sigma = nu_sigma2
      ))
      # Accept draws
      accept <- ll1 - ll0
      accept <- ifelse(test = accept > 0, yes = 0, no = accept)
      accept <- ifelse(test = runif(
        nrow(nu0)) < exp(x = accept),
        yes = 1,
        no = 0
        )
      nu1[which(x = accept != 1), ] <- nu0[which(x = accept != 1), ]
      if (iter_num < tune_int) {
        recent_nu_draws <- abind::abind(
          recent_nu_draws,
          array(data = t(nu1), dim = c(1, 1, nrow(nu0))),
          along = 1
        )
      } else {
        recent_nu_draws <- abind::abind(
          array(
            data = recent_nu_draws[-1, , ],
            dim = dim(x = recent_nu_draws) - c(1, 0, 0)
          ),
          array(data = t(nu1), dim = c(1, 1, nrow(nu0))),
          along = 1
        )
      }
      # Adjust tuning rate and candidate covariance matrix
      nurate <- apply(
        X = recent_nu_draws,
        MARGIN = c(3),
        FUN = function(x) {
          mean(x = !duplicated(x))
        }
      )
      if (iter_num %in% seq(from = min_tune, to = max_tune, by = tune_int)) {
        n.scale <- (n.scale * qnorm(p = .23 / 2)) /
          qnorm(p = (nurate * .997 + .001) / 2)
        cand.n.var <- (
          ((iter_num / (max_tune + weight * max_tune)) * array(
            data = apply(
              X = array(
                data = recent_nu_draws[
                  seq(from = 1,
                      to = tune_int,
                      by = ceiling(x = tune_int / 10)
                  ), , ],
                dim = c(
                  length(
                    x = seq(
                      from = 1,
                      to = tune_int,
                      by = ceiling(x = tune_int / 10)
                    )
                  ),
                  ncol(nu0), nrow(nu0)
                )
              ),
              MARGIN = c(3),
              FUN = var),
            dim = c(ncol(nu0), ncol(nu0), nrow(nu0))
          )) +
            ((1 - (iter_num / (max_tune + weight * max_tune)))) * array(
              data = unlist(
                lapply(
                  X = n.scale,
                  function(x) {
                    diag(x = x, ncol(nu0))
                  }
                )
              ),
              dim = c(ncol(nu0), ncol(nu0), nrow(nu0)))
        )
      }
      if (iter_num > burn & iter_num %in% seq(
        from = burn + 1,
        to = niter,
        by = thin)
      ) {
        nu_draws <- abind::abind(
          nu_draws,
          array(data = t(nu1), dim = c(1, 1, nrow(nu0))),
          along = 1
        )
      }
    } else {
      nu1 <- nu0
    }
    # Draw zetas
    if (est_zeta) {
      zeta1 <- zeta0 + matrix(
        data = sapply(
          X = seq_len(nrow(zeta0)),
          FUN = function(i) {
            MASS::mvrnorm(
              n = 1,
              mu = rep(0, ncol(zeta0)),
              Sigma = cand_z_var[, , i])
          }
        ),
        nrow = nrow(zeta0),
        ncol = ncol(zeta0),
        byrow = T
      )
      p0 <- obj_fun(y = y, nu = array(data = nu0, dim = dim(y)),
                    lambda = lambda, gamma = gamma, omega = omega0,
                    zeta = zeta0)$p
      tmp0 <- rowSums(x = log((p0^y) * (1 - p0) ^ (1 - y)), na.rm = T)
      ll0 <- tmp0 + log(x = mvtnorm::dmvnorm(
        x = zeta0,
        mean = zeta_mu,
        sigma = zeta_sigma2
      ))
      p1 <- obj_fun(y = y, nu = array(data = nu0, dim = dim(y)),
                    lambda = lambda, gamma = gamma, omega = omega0,
                    zeta = zeta1)$p
      tmp1 <- rowSums(x = log((p1^y) * (1 - p1) ^ (1 - y)), na.rm = T)
      ll1 <- tmp1 + log(x = mvtnorm::dmvnorm(
        x = zeta1,
        mean = zeta_mu,
        sigma = zeta_sigma2
      ))
      # Accept draws
      accept <- ll1 - ll0
      accept <- ifelse(test = accept > 0, yes = 0, no = accept)
      accept <- ifelse(
        test = runif(nrow(zeta0)) < exp(x = accept),
        yes = 1,
        no = 0
      )
      zeta1[which(x = accept != 1), ] <- zeta0[which(x = accept != 1), ]
      if (iter_num < tune_int) {
        recent_zeta_draws <- abind::abind(
          recent_zeta_draws,
          array(data = t(zeta1), c(1, ncol(zeta0), nrow(zeta0))),
          along = 1
        )
      } else {
        recent_zeta_draws <- abind::abind(
          array(
            recent_zeta_draws[-1, , ],
            dim(x = recent_zeta_draws) - c(1, 0, 0)
          ),
          array(t(x = zeta1), c(1, ncol(zeta0), nrow(zeta0))),
          along = 1
        )
      }
      # Adjust tuning rate and candidate covariance matrix
      zetarate <- apply(
        X = recent_zeta_draws,
        MARGIN = c(3),
        FUN = function(x) {
          mean(x = !duplicated(x = x))
        }
      )
      if (iter_num %in% seq(from = min_tune, to = max_tune, by = tune_int)) {
        z.scale <- (z.scale * qnorm(p = .23 / 2)) /
          qnorm(p = (zetarate * .997 + .001) / 2)
        cand_z_var <- (
          ((iter_num / (max_tune + weight * max_tune)) * array(
            data = apply(
              X = array(
                data = recent_zeta_draws[seq(
                  from = 1,
                  to = tune_int,
                  by = ceiling(x = tune_int / 10)
                ), , ],
                dim = c(length(x = seq(
                  from = 1,
                  to = tune_int,
                  by = ceiling(x = tune_int / 10))
                ), ncol(zeta0), nrow(zeta0))
              ),
              MARGIN = c(3),
              FUN = var),
            dim = c(ncol(zeta0), ncol(zeta0), nrow(zeta0))
          )) +
            ((1 - (iter_num / (max_tune + weight * max_tune)))) * array(
              data = unlist(lapply(
                X = z.scale,
                FUN = function(x) {
                  diag(x = x, ncol(zeta0))
                }
              )),
              dim = c(ncol(zeta0), ncol(zeta0), nrow(zeta0)))
        )
      }
      if (iter_num > burn & iter_num %in% seq(
        from = burn + 1,
        to = niter,
        by = thin)
      ) {
        zeta_draws <- abind::abind(
          zeta_draws, array(
            data = t(zeta1),
            dim = c(1, ncol(zeta0), nrow(zeta0))
          ),
          along = 1
        )
      }
    } else {
      zeta1 <- zeta0
    }
    #Update logLikelihood and estimates
    ll <- obj_fun(y = y, nu = array(data = nu1, dim = dim(y)), lambda = lambda,
                 gamma = gamma, omega = omega1, zeta = zeta1)$loglikelihood
    omega0 <- omega1
    nu0 <- nu1
    zeta0 <- zeta1
    if (i %in% seq(from = 1, to = niter, by = thin) & verbose_mcmh) {
      cat("Iter.",
          format(x = i, width = 4, justify = "left"), sep = "",
          append = TRUE
      )
      cat(" ll =",
          format(x = round(x = ll, digits = 2), nsmall = 2, width = 10),
          sep = "",
          append = TRUE
      )
      if (est_omega) {
        cat(" Mean omega Rate = ",
            format(x = round(x = 100 * mean(x = omegarate)), width = 3),
            "% ",
            "Mean o_scale = ",
            format(x = round(x = mean(x = o_scale), 3), width = 6, nsmall = 3),
            sep = "",
            append = TRUE
        )
      }
      if (est_nu) {
        cat(" Mean nu Rate = ",
            format(x = round(x = 100 * mean(x = nurate)), width = 3),
            "% ",
            "Mean n.scale = ",
            format(x = round(x = mean(x = n.scale), digits = 3),
                   width = 6,
                   nsmall = 3
            ),
            sep = "",
            append = TRUE
        )
      }
      if (est_zeta) {
        cat(" Mean zeta Rate = ",
            format(x = round(x = 100 * mean(x = zetarate)), width = 3),
            "% ",
            "Mean z.scale = ",
            format(x = round(x = mean(x = z.scale), digits = 3),
                   width = 6,
                   nsmall = 3
            ),
            sep = "",
            append = TRUE)
      }
      cat(" Time = ",
          format(x = Sys.time(), format = "%H:%M:%S"),
          sep = "",
          append = TRUE
      )
      cat("\n", sep = ",", append = TRUE)
    }
  }
  return(list(
    "omega_draws" = omega_draws,
    "nu_draws" = nu_draws,
    "zeta_draws" = zeta_draws,
    "cand_o_var" = cand_o_var,
    "cand_n_var" = cand.n.var,
    "cand_z_var" = cand_z_var
  ))
}
