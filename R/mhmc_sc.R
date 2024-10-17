#-------------------------------------------------------------------------------
#' MHMC Parameter Estimates for Single Chain
#'
#' This function uses the Metropolis-Hastings algorithm for a single chain to
#' calculate parameter estimates using the Markov Chain Monte Carlo method. The
#' method implemented follows Patz and Junker (1999). The approach to tuning the
#' scale and covariance matrix follows BDA and SAS 9.2 User Guide,
#' 2nd Ed. "The MCMC Procedure: Tuning the Proposal Distribution".
#'
#' @param y Matrix of item responses (K by IJ).
#' @param obj_fun A function that calculates predictions and log-likelihood
#' values for the selected model (character).
#' @param link Choose between "logit" or "probit" link functions.
#' @param est_omega Determines whether omega is estimated (logical).
#' @param est_lambda Determines whether nu is estimated (logical).
#' @param est_zeta Determines whether zeta is estimated (logical).
#' @param est_nu Determines whether nu is estimated (logical).
#' @param omega0 Starting or known values for omega (K by MN).
#' @param gamma0 Starting or known values for gamma (JM by MN).
#' @param lambda0 Starting or known values for lambda (IJ by JM).
#' @param zeta0 Starting or known values for  zeta (K by JM).
#' @param nu0 Starting or known values for nu (IJ by 1).
#' @param kappa0 Starting or known values for kappa (K by IJ).
#' @param omega_mu Mean prior for omega (1 by MN).
#' @param omega_sigma2 Covariance prior for omega (MN by MN).
#' @param lambda_mu Mean prior for lambda (1 by JM)
#' @param lambda_sigma2 Covariance prior for lambda (JM by JM)
#' @param zeta_mu Mean prior for zeta (1 by JM).
#' @param zeta_sigma2 Covariance prior for zeta (JM by JM).
#' @param nu_mu Mean prior for nu (1 by 1).
#' @param nu_sigma2 Covariance prior for nu (1 by 1).
#' @param burn Number of iterations at the beginning of an MCMC run to discard
#' (scalar).
#' @param thin Determines every nth observation retained (scalar).
#' @param min_tune Determines when tunning begins (scalar).
#' @param tune_int MHMC tuning interval (scalar).
#' @param max_tune Determines when tunning ends (scalar).
#' @param niter Number of iterations of the MHMC sampler.
#' @param weight Determines the weight of old versus new covariance matrix.
#' @param verbose_mhmc Print progress of MHMC sampler.
#'
#' @return List with elements omega_draws (list of (niter - burn) / thin draws
#' for K by MN omega matrix), lambda_draws (list of (niter - burn) / thin draws
#' for IJ by JM lambda matrix), zeta_draws (list of (niter - burn) / thin draws
#' for K by JM zeta matrix), nu_draws (list of (niter - burn) / thin draws
#' for IJ by 1 nu matrix), cand_o_var (list of K final MN by MN candidate
#' proposal covariance matrices for omega for each examinee), cand_l_var (list
#' of IJ final JM by JM candidate proposal covariance matrices for lambda for
#' each item), cand_z_var (list of final JM by JM candidate proposal covariance
#' matrices for zeta for all examinees), and cand_n_var (list of IJ final scalar
#' candidate proposal variances for nu for all items).
#'
#' @references
#'
#' Patz, R. J., & Junker, B. W. (1999). A Straightforward Approach to Markov
#' Chain Monte Carlo Methods for Item Response Models. Journal of Educational
#' and Behavioral Statistics, 24(2), 146.
#'
#' @keywords internal
#-------------------------------------------------------------------------------

mhmc_sc <- function(
    y = y,
    obj_fun = NULL,
    link = NULL,
    est_omega = TRUE,
    est_lambda = TRUE,
    est_zeta = TRUE,
    est_nu = TRUE,
    omega0 = NULL,
    gamma0 = NULL,
    lambda0 = NULL,
    zeta0 = NULL,
    nu0 = NULL,
    kappa0 = NULL,
    omega_mu = NULL,
    omega_sigma2 = NULL,
    lambda_mu = NULL,
    lambda_sigma2 = NULL,
    zeta_mu = NULL,
    zeta_sigma2 = NULL,
    nu_mu = NULL,
    nu_sigma2 = NULL,
    burn = NULL,
    thin = NULL,
    min_tune = NULL,
    tune_int = NULL,
    max_tune = NULL,
    niter = NULL,
    weight = 1,
    verbose_mhmc = FALSE
) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package \"MASS\" needed for the mhmc_sc function to work. Please
         install.",
         call. = FALSE)
  }
  if (!requireNamespace("abind", quietly = TRUE)) {
    stop("Package \"abind\" needed for the mhmc_sc function to work. Please
         install.",
         call. = FALSE)
  }
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package \"mvtnorm\" needed for the mhmc_sc function to work. Please
         install.",
         call. = FALSE)
  }

  if (verbose_mhmc) {
    cat(
      "\n",
      "MCMC LOG",
      "Start Time",
      format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
      "\n\n",
      sep = " "
    )
  }
  set.seed(1)
  # Initialize objects to save recent draws for tuning, objects to save final
  # draws, and objects to save final draws. For candidate proposal variance(see
  # BDA3 p. 296.
  if (est_omega) {
    recent_omega_draws <- array(
      data = t(omega0),
      dim = c(1, ncol(x = omega0), nrow(x = omega0))
    )
    omega_draws <- NULL
    o_scale <- 2.38 / sqrt(x = ncol(x = omega0))
    o_prop <- array(
      data = diag(x = .75, nrow = ncol(x = omega0)),
      dim = c(ncol(x = omega0), ncol(x = omega0), nrow(x = omega0))
    )
    cand_o_var <- o_prop * o_scale
  } else {
    omega_draws <- NULL
    cand_o_var <- NULL
  }
  if (est_lambda) {
    recent_lambda_draws <- array(
      data = t(lambda0),
      dim = c(1, ncol(x = lambda0), nrow(x = lambda0))
    )
    lambda_draws <- NULL
    l_scale <- 2.38 / sqrt(x = ncol(x = lambda0))
    l_prop <- array(
      data = diag(x = .75, nrow = ncol(lambda0)),
      dim = c(ncol(x = lambda0), ncol(x = lambda0), nrow(x = lambda0))
    )
    cand_l_var <- l_prop * l_scale
  } else {
    lambda_draws <- NULL
    cand_l_var <- NULL
  }
  if (est_nu) {
    recent_nu_draws <- array(
      data = t(nu0),
      dim = c(1, 1, nrow(x = nu0))
    )
    nu_draws <- NULL
    n_scale <- 2.38 / sqrt(x = ncol(x = nu0))
    n_prop <- array(
      data = diag(x = .25, nrow = 1),
      dim = c(1, 1, nrow(x = nu0))
    )
    cand_n_var <- n_prop * n_scale
  } else {
    nu_draws <- NULL
    cand_n_var <- NULL
  }
  if (est_zeta) {
    recent_zeta_draws <- array(
      data = t(zeta0),
      dim = c(1, ncol(x = zeta0), nrow(x = zeta0))
    )
    zeta_draws <- NULL
    z_scale <- 2.38 / sqrt(x = ncol(x = zeta0))
    z_prop <- array(
      data = diag(x = .05, nrow = ncol(x = zeta0)),
      dim = c(ncol(x = zeta0), ncol(x = zeta0), nrow(x = zeta0))
    )
    cand_z_var <- z_prop * z_scale
  } else {
    zeta_draws <- NULL
    cand_z_var <- NULL
  }
  for (i in 1:niter) {
    iter_num <- i
    # Draw omegas
    if (est_omega) {
      omega1 <- omega0 + matrix(
        data = sapply(
          X = seq_len(length.out = nrow(x = omega0)),
          FUN = function(i) {
            MASS::mvrnorm(
              n = 1,
              mu = rep(0, ncol(x = omega0)),
              Sigma = cand_o_var[, , i])
          }
        ),
        nrow = nrow(x = omega0),
        ncol = ncol(x = omega0),
        byrow = TRUE
      )
      p0 <- obj_fun(y = y, omega = omega0, gamma = gamma0, lambda = lambda0,
                    zeta = zeta0, nu = nu0, kappa = kappa0, link = link)$p
      tmp0 <- rowSums(x = log((p0^y) * (1 - p0) ^ (1 - y)), na.rm = TRUE)
      ll0 <- tmp0 + log(x = mvtnorm::dmvnorm(
        x = omega0,
        mean = omega_mu,
        sigma = omega_sigma2
      ))
      p1 <- obj_fun(y = y, omega = omega1, gamma = gamma0, lambda = lambda0,
                    zeta = zeta0, nu = nu0, kappa = kappa0, link = link)$p
      tmp1 <- rowSums(x = log((p1^y) * (1 - p1) ^ (1 - y)), na.rm = TRUE)
      ll1 <- tmp1 + log(x = mvtnorm::dmvnorm(
        x = omega1,
        mean = omega_mu,
        sigma = omega_sigma2
      ))
      # Accept draws
      accept <- ll1 - ll0
      accept <- ifelse(test = accept > 0, yes = 0, no = accept)
      accept <- ifelse(
        test = runif(n = nrow(x = omega0)) < exp(x = accept),
        yes = 1,
        no = 0
      )
      omega1[which(x = accept != 1), ] <- omega0[which(x = accept != 1), ]
      if (iter_num < tune_int) {
        recent_omega_draws <- abind::abind(
          recent_omega_draws,
          array(
            data = t(omega1),
            dim = c(1, ncol(x = omega0), nrow(x = omega0))
          ),
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
                      ncol(x = omega0),
                      nrow(x = omega0)
                    )
                  ),
                  MARGIN = c(3),
                  FUN = var
                ),
                dim = c(ncol(x = omega0), ncol(x = omega0), nrow(x = omega0))
              )
          ) +
            ((1 - (iter_num / (max_tune + weight * max_tune)))) *
            array(
              data = unlist(lapply(X = o_scale, function(x) {
                diag(x = x, ncol(x = omega0))
              })),
              dim = c(ncol(x = omega0), ncol(x = omega0), nrow(x = omega0))
            )
        )
      }
      if (iter_num > burn && iter_num %in% seq(
        from = burn + 1,
        to = niter,
        by = thin
      )
      ) {
        omega_draws <- abind::abind(
          omega_draws,
          array(
            data = t(omega1),
            dim = c(1, ncol(x = omega0), nrow(x = omega0))
          ),
          along = 1)
      }
    } else {
      omega1 <- omega0
    }
    # Draw lambdas
    if (est_lambda) {
      lambda1 <- lambda0 + matrix(
        data = sapply(
          X = seq_len(length.out = nrow(x = lambda0)),
          FUN = function(i) {
            MASS::mvrnorm(
              n = 1,
              mu = rep(0, ncol(x = lambda0)),
              Sigma = cand_l_var[, , i])
          }
        ),
        nrow = nrow(x = lambda0),
        ncol = ncol(x = lambda0),
        byrow = TRUE
      )
      p0 <- obj_fun(y = y, omega = omega0, gamma = gamma0, lambda = lambda0,
                    zeta = zeta0, nu = nu0, kappa = kappa0, link = link)$p
      tmp0 <- colSums(x = log((p0^y) * (1 - p0) ^ (1 - y)), na.rm = TRUE)
      ll0 <- tmp0 + log(x = mvtnorm::dmvnorm(
        x = lambda0,
        mean = lambda_mu,
        sigma = lambda_sigma2
      ))
      p1 <- obj_fun(y = y, omega = omega0, gamma = gamma0, lambda = lambda1,
                    zeta = zeta0, nu = nu0, kappa = kappa0, link = link)$p
      tmp1 <- colSums(x = log(x = (p1^y) * (1 - p1) ^ (1 - y)), na.rm = TRUE)
      ll1 <- tmp1 + log(x = mvtnorm::dmvnorm(
        x = lambda1,
        mean = lambda_mu,
        sigma = lambda_sigma2
      ))
      # Accept draws
      accept <- ll1 - ll0
      accept <- ifelse(test = accept > 0, yes = 0, no = accept)
      accept <- ifelse(
        test = runif(nrow(x = lambda0)) < exp(x = accept),
        yes = 1,
        no = 0
      )
      lambda1[which(x = accept != 1), ] <- lambda0[which(x = accept != 1), ]
      if (iter_num < tune_int) {
        recent_lambda_draws <- abind::abind(
          recent_lambda_draws,
          array(
            data = t(lambda1),
            dim = c(1, ncol(x = lambda0), nrow(x = lambda0))
          ),
          along = 1
        )
      } else {
        recent_lambda_draws <- abind::abind(
          array(
            recent_lambda_draws[-1, , ],
            dim(x = recent_lambda_draws) - c(1, 0, 0)
          ),
          array(
            data = t(x = lambda1),
            dim = c(1, ncol(x = lambda0), nrow(x = lambda0))
          ),
          along = 1
        )
      }
      # Adjust tuning rate and candidate covariance matrix
      lambdarate <- apply(
        X = recent_lambda_draws,
        MARGIN = c(3),
        FUN = function(x) {
          mean(x = !duplicated(x = x))
        }
      )
      if (iter_num %in% seq(from = min_tune, to = max_tune, by = tune_int)) {
        l_scale <- (l_scale * qnorm(p = .23 / 2)) /
          qnorm(p = (lambdarate * .997 + .001) / 2)
        cand_l_var <- (
          (
            (iter_num / (max_tune + weight * max_tune)) *
              array(
                data = apply(
                  X = array(
                    data = recent_lambda_draws[
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
                      ncol(x = lambda0),
                      nrow(x = lambda0)
                    )
                  ),
                  MARGIN = c(3),
                  FUN = var
                ),
                dim = c(ncol(x = lambda0), ncol(x = lambda0), nrow(x = lambda0))
              )
          ) +
            ((1 - (iter_num / (max_tune + weight * max_tune)))) *
            array(
              data = unlist(lapply(X = l_scale, function(x) {
                diag(x = x, ncol(x = lambda0))
              })),
              dim = c(ncol(x = lambda0), ncol(x = lambda0), nrow(x = lambda0))
            )
        )
      }
      if (iter_num > burn && iter_num %in% seq(
        from = burn + 1,
        to = niter,
        by = thin
      )
      ) {
        lambda_draws <- abind::abind(
          lambda_draws,
          array(
            data = t(lambda1),
            dim = c(1, ncol(x = lambda0), nrow(x = lambda0))
          ),
          along = 1)
      }
    } else {
      lambda1 <- lambda0
    }
    # Draw zetas
    if (est_zeta) {
      zeta1 <- zeta0 + matrix(
        data = sapply(
          X = seq_len(length.out = nrow(x = zeta0)),
          FUN = function(i) {
            MASS::mvrnorm(
              n = 1,
              mu = rep(0, ncol(x = zeta0)),
              Sigma = cand_z_var[, , i])
          }
        ),
        nrow = nrow(x = zeta0),
        ncol = ncol(x = zeta0),
        byrow = TRUE
      )
      p0 <- obj_fun(y = y, omega = omega0, gamma = gamma0, lambda = lambda0,
                    zeta = zeta0, nu = nu0, kappa = kappa0, link = link)$p
      tmp0 <- rowSums(x = log((p0^y) * (1 - p0) ^ (1 - y)), na.rm = TRUE)
      ll0 <- tmp0 + log(x = mvtnorm::dmvnorm(
        x = zeta0,
        mean = zeta_mu,
        sigma = zeta_sigma2
      ))
      p1 <- obj_fun(y = y, omega = omega0, gamma = gamma0, lambda = lambda0,
                    zeta = zeta1, nu = nu0, kappa = kappa0, link = link)$p
      tmp1 <- rowSums(x = log((p1^y) * (1 - p1) ^ (1 - y)), na.rm = TRUE)
      ll1 <- tmp1 + log(x = mvtnorm::dmvnorm(
        x = zeta1,
        mean = zeta_mu,
        sigma = zeta_sigma2
      ))
      # Accept draws
      accept <- ll1 - ll0
      accept <- ifelse(test = accept > 0, yes = 0, no = accept)
      accept <- ifelse(
        test = runif(nrow(x = zeta0)) < exp(x = accept),
        yes = 1,
        no = 0
      )
      zeta1[which(x = accept != 1), ] <- zeta0[which(x = accept != 1), ]
      if (iter_num < tune_int) {
        recent_zeta_draws <- abind::abind(
          recent_zeta_draws,
          array(data = t(zeta1), c(1, ncol(x = zeta0), nrow(x = zeta0))),
          along = 1
        )
      } else {
        recent_zeta_draws <- abind::abind(
          array(
            recent_zeta_draws[-1, , ],
            dim(x = recent_zeta_draws) - c(1, 0, 0)
          ),
          array(t(x = zeta1), c(1, ncol(x = zeta0), nrow(x = zeta0))),
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
        z_scale <- (z_scale * qnorm(p = .23 / 2)) /
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
                ), ncol(x = zeta0), nrow(x = zeta0))
              ),
              MARGIN = c(3),
              FUN = var),
            dim = c(ncol(x = zeta0), ncol(x = zeta0), nrow(x = zeta0))
          )) +
            ((1 - (iter_num / (max_tune + weight * max_tune)))) * array(
              data = unlist(lapply(
                X = z_scale,
                FUN = function(x) {
                  diag(x = x, ncol(x = zeta0))
                }
              )),
              dim = c(ncol(x = zeta0), ncol(x = zeta0), nrow(x = zeta0)))
        )
      }
      if (iter_num > burn && iter_num %in% seq(
        from = burn + 1,
        to = niter,
        by = thin
      )
      ) {
        zeta_draws <- abind::abind(
          zeta_draws, array(
            data = t(zeta1),
            dim = c(1, ncol(x = zeta0), nrow(x = zeta0))
          ),
          along = 1
        )
      }
    } else {
      zeta1 <- zeta0
    }
    # Draw nus
    if (est_nu) {
      nu1 <- nu0 + matrix(
        data = sapply(
          X = seq_len(length.out = nrow(x = nu0)),
          FUN = function(j) {
            MASS::mvrnorm(
              n = 1,
              mu = rep(0, 1),
              Sigma = cand_n_var[, , j]
            )
          }
        ), nrow = nrow(x = nu0),
        ncol = ncol(x = nu0),
        byrow = TRUE
      )
      p0 <- obj_fun(y = y, omega = omega0, gamma = gamma0, lambda = lambda0,
                    zeta = zeta0, nu = nu0, kappa = kappa0, link = link)$p
      tmp0 <- colSums(x = log((p0^y) * (1 - p0) ^ (1 - y)), na.rm = TRUE)
      ll0 <- tmp0 + log(x = mvtnorm::dmvnorm(
        x = nu0,
        mean = nu_mu,
        sigma = nu_sigma2
      ))
      p1 <- obj_fun(y = y, omega = omega0, gamma = gamma0, lambda = lambda0,
                    zeta = zeta0, nu = nu1, kappa = kappa0, link = link)$p
      tmp1 <- colSums(x = log((p1 ^ y) * (1 - p1) ^ (1 - y)), na.rm = TRUE)
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
          array(data = t(nu1), dim = c(1, 1, nrow(x = nu0))),
          along = 1
        )
      } else {
        recent_nu_draws <- abind::abind(
          array(
            data = recent_nu_draws[-1, , ],
            dim = dim(x = recent_nu_draws) - c(1, 0, 0)
          ),
          array(data = t(nu1), dim = c(1, 1, nrow(x = nu0))),
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
        n_scale <- (n_scale * qnorm(p = .23 / 2)) /
          qnorm(p = (nurate * .997 + .001) / 2)
        cand_n_var <- (
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
                  ncol(x = nu0), nrow(x = nu0)
                )
              ),
              MARGIN = c(3),
              FUN = var),
            dim = c(ncol(x = nu0), ncol(x = nu0), nrow(x = nu0))
          )) +
            ((1 - (iter_num / (max_tune + weight * max_tune)))) * array(
              data = unlist(
                lapply(
                  X = n_scale,
                  function(x) {
                    diag(x = x, ncol(x = nu0))
                  }
                )
              ),
              dim = c(ncol(x = nu0), ncol(x = nu0), nrow(x = nu0)))
        )
      }
      if (iter_num > burn && iter_num %in% seq(
        from = burn + 1,
        to = niter,
        by = thin
      )
      ) {
        nu_draws <- abind::abind(
          nu_draws,
          array(data = t(nu1), dim = c(1, 1, nrow(x = nu0))),
          along = 1
        )
      }
    } else {
      nu1 <- nu0
    }
    #Update logLikelihood and estimates
    ll <- obj_fun(y = y, omega = omega1, gamma = gamma0, lambda = lambda1,
                  zeta = zeta1, nu = nu1, kappa = kappa0,
                  link = link)$loglikelihood
    omega0 <- omega1
    lambda0 <- lambda1
    zeta0 <- zeta1
    nu0 <- nu1
    if (i %in% seq(from = 1, to = niter, by = thin) && verbose_mhmc) {
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
        cat(" Mean omega rate =",
            format(x = round(x = 100 * mean(x = omegarate)), width = 3),
            "% ",
            "Mean o_scale =",
            format(x = round(x = mean(x = o_scale), 2),
                   width = 5,
                   nsmall = 2
            ),
            sep = "",
            append = TRUE
        )
      }
      if (est_lambda) {
        cat(" Mean lambda rate =",
            format(x = round(x = 100 * mean(x = omegarate)), width = 3),
            "% ",
            "Mean l_scale =",
            format(x = round(x = mean(x = o_scale), 2),
                   width = 5,
                   nsmall = 2
            ),
            sep = "",
            append = TRUE
        )
      }
      if (est_zeta) {
        cat(" Mean zeta rate =",
            format(x = round(x = 100 * mean(x = zetarate)), width = 3),
            "% ",
            "Mean z_scale =",
            format(x = round(x = mean(x = z_scale), digits = 2),
                   width = 5,
                   nsmall = 2
            ),
            sep = "",
            append = TRUE)
      }
      if (est_nu) {
        cat(" Mean nu rate =",
            format(x = round(x = 100 * mean(x = nurate)), width = 3),
            "% ",
            "Mean n_scale =",
            format(x = round(x = mean(x = n_scale), digits = 2),
                   width = 5,
                   nsmall = 2
            ),
            sep = "",
            append = TRUE
        )
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
    "lambda_draws" = lambda_draws,
    "zeta_draws" = zeta_draws,
    "nu_draws" = nu_draws,
    "cand_o_var" = cand_o_var,
    "cand_l_var" = cand_l_var,
    "cand_z_var" = cand_z_var,
    "cand_n_var" = cand_n_var
  ))
}
