#-------------------------------------------------------------------------------
#' MHRM Parameter Estimates for Multiple Chains
#'
#' This function calculates mhrm parameter estimates for multiple chains.
#'
#' @param y Matrix of item responses (K by IJ).
#' @param obj_fun A function that calculates predictions and log-likelihood
#' values for the selected model (character).
#' @param link Choose between "logit" or "probit" link functions.
#' @param est_omega Determines whether omega is estimated (logical).
#' @param est_lambda Determines whether lambda is estimated (logical).
#' @param est_nu Determines whether nu is estimated (logical).
#' @param est_zeta Determines whether zeta is estimated (logical).
#' @param lambda0 Matrix of item structure parameters (IJ by JM).
#' @param kappa0 Matrix of item guessing parameters (IJ by 1).
#' @param gamma0 Either a matrix of experimental structure parameters (JM by MN)
#' or the name in quotes of the desired R stats contrast function (i.e.,
#' "contr.helmert", "contr.poly", "contr.sum", "contr.treatment", or
#' "contr.SAS"). If using the R stats contrast function the user must also
#' specify J, M, and N, as well as ensure that items in y are arranged so that
#' the first set of I items correspond to the first level if the contrast, the
#' next set of I items correspond to the second level of the contrast, etc. For
#' example, in an experiment with two conditions (i.e., J = 2) where the user
#' requests two contrasts (i.e., N = 2) from the "contr.treatment" function, the
#' first set of I items will all receive a contrast code of 0 and the second set
#' of I items will all receive a contrast code of 1. In an experiment with three
#' conditions (i.e., J = 3) where the user requests three contrasts (i.e.,
#' N = 3) from the "contr.poly" function, first set of I items will receive the
#' lowest value code for linear and quadratic contrasts, the second set of I
#' items will all receive the middle value code for linear and quadratic
#' contrasts, and the last set of I items will all receive the highest value
#' code for linear and quadratic contrasts.
#' @param omega0 Starting values for omega.
#' @param nu0 Starting values for nu (IJ by 1).
#' @param zeta0 Starting values for zeta.
#' @param omega_mu Vector of means prior for omega (1 by MN).
#' @param omega_sigma2 Covariance matrix prior for omega (MN by MN).
#' @param lambda_mu Mean prior for lambda (1 by JM)
#' @param lambda_sigma2 Covariance prior for lambda (JM by JM)
#' @param zeta_mu Vector of means prior for zeta (1 by JM).
#' @param zeta_sigma2 Covariance matrix prior for zeta (JM by JM).
#' @param nu_mu Prior mean for nu (scalar).
#' @param nu_sigma2 Prior variance for nu (scalar).
#' @param constraints Item parameter constraints.
#' @param J Number of conditions (required if using the R stats contrast
#' function).
#' @param M Number of ability (or trait) dimensions (required if using the R
#' stats contrast function).
#' @param N Number of contrasts including intercept (required if using the R
#' stats contrast function).
#' @param ... Additional arguments.
#'
#' @return List with elements for all parameters estimated, information values
#' for all parameters estimated, and the model log-likelihood value.
#'
#' @references
#' Cai, L. (2010). High-dimensional exploratory item factor analysis by a
#' Metropolis-Hastings Robbins-Monro algorithm. \emph{Psychometrika, 75(1)},
#' 33-57.
#'
#' Cai, L. (2010). Metropolis-Hastings Robbins-Monro algorithm for confirmatory
#' item factor analysis. \emph{Journal of Educational and Behavioral Statistics,
#' 35(3)}, 307-335.
#'
#' @examples
#'# mhrmfit <- mhrm(y = ex3$y, obj_fun = dich_response_model, est_omega = TRUE,
#'#                 est_lambda = FALSE, est_zeta = TRUE, est_nu = TRUE,
#'#                 omega0 = array(data = 0, dim = dim(ex3$omega)),
#'#                 gamma0 = ex3$gamma, lambda0 = ex3$lambda,
#'#                 zeta0 = array(data = 0, dim = dim(ex3$zeta)),
#'#                 nu0 = array(data = 0, dim = dim(ex3$nu)),
#'#                 kappa0 = ex3$kappa, omega_mu = ex3$omega_mu,
#'#                 omega_sigma2 = ex3$omega_sigma2, zeta_mu = ex3$zeta_mu,
#'#                 zeta_sigma2 = ex3$zeta_sigma2, nu_mu = ex3$nu_mu,
#'#                 nu_sigma2 = ex3$nu_sigma2)
#'# summary(mhrmfit)
#'# plot(mhrmfit)
#'
#' @keywords internal
#-------------------------------------------------------------------------------

mhrm <- function(
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
    constraints = NULL,
    J = NULL,
    M = NULL,
    N = NULL,
    ...
) {
  if (!requireNamespace("abind", quietly = TRUE)) {
    stop("Package \"abind\" needed for the mhrm function to work. Please
         install.",
         call. = FALSE)
  }
  ellipsis <- list(...)
  if (is.null(x = ellipsis$chains)) {
    chains <- 3
  } else {
    chains <- ellipsis$chains
  }
  if (is.null(x = ellipsis$burn)) {
    burn <- 0
  } else {
    burn <- ellipsis$burn
  }
  if (is.null(x = ellipsis$thin)) {
    thin <- 5
  } else {
    thin <- ellipsis$thin
  }
  if (is.null(x = ellipsis$min_tune)) {
    min_tune <- 0
  } else {
    min_tune <- ellipsis$min_tune
  }
  if (is.null(x = ellipsis$tune_int)) {
    tune_int <- 0
  } else {
    tune_int <- ellipsis$tune_int
  }
  if (is.null(x = ellipsis$max_tune)) {
    max_tune <- 0
  } else {
    max_tune <- ellipsis$max_tune
  }
  if (is.null(x = ellipsis$niter)) {
    niter <- 6
  } else {
    niter <- ellipsis$niter
  }
  if (is.null(x = ellipsis$max_iter_mhrm)) {
    max_iter_mhrm <- 200
  } else {
    max_iter_mhrm <- ellipsis$max_iter_mhrm
  }
  if (is.null(x = ellipsis$verbose_mhrm)) {
    verbose_mhrm <- TRUE
  } else {
    verbose_mhrm <- ellipsis$verbose_mhrm
  }
  if (is.character(x = gamma0)) {
    if (gamma0 %in% c("contr.helmert", "contr.poly", "contr.sum",
                      "contr.treatment", "contr.SAS")) {
      contrast_codes <- cbind(1, get(x = gamma0)(n = J))[, 1:N]
      tmp <- matrix(data = 0, nrow = J * M, ncol = M * N)
      for (j in 1:J) {
        for (m in 1:M) {
          tmp[(m + M * (j - 1)), (((m - 1) * N + 1):((m - 1) * N + N))] <-
            contrast_codes[j, ]
        }
      }
      gamma0 <- tmp
    }
  }
  # STEP 0: MCMC burn in -------------------------------------------------------
  if (verbose_mhrm) {
    cat(
      "MCMC Burn-In Start Time",
      format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
      "\n",
      sep = " "
    )
  }
  mhmcburn <- mhmc_mc(
    chains = chains, y = y, obj_fun = obj_fun, link = link,
    est_omega = est_omega, est_lambda = est_lambda, est_zeta = est_zeta,
    est_nu = est_nu, omega0 = omega0, gamma0 = gamma0, lambda0 = lambda0,
    zeta0 = zeta0, nu0 = nu0, kappa0 = kappa0, omega_mu = omega_mu,
    omega_sigma2 = omega_sigma2, lambda_mu = lambda_mu,
    lambda_sigma2 = lambda_sigma2, zeta_mu = zeta_mu, zeta_sigma2 = zeta_sigma2,
    nu_mu = nu_mu, nu_sigma2 = nu_sigma2, burn = burn,
    thin = thin, min_tune = min_tune, tune_int = tune_int, max_tune = max_tune,
    niter = niter
  )

  # Update initial estimates and variance of candidates
  if (est_omega) {
    omega0 <- mhmcburn$omegaEAP
    info0_omega <- diag(x = 10, nrow = ncol(omega0))
  }
  if (est_lambda) {
    lambda0 <- mhmcburn$lambdaEAP
    if (!is.null(constraints)) {
      lambda0 <- array(data = t(c(lambda0)) %*% constraints[[4]],
                       dim = dim(lambda0))
    }
    info0_lambda <- diag(x = 10, nrow = ncol(lambda0))
  }
  if (est_zeta) {
    zeta0 <- mhmcburn$zetaEAP
  }
  if (est_nu) {
    nu0 <- mhmcburn$nuEAP
    if (!is.null(constraints)) {
      nu0 <- t(t(nu0) %*% constraints[[2]])
    }
    info0_nu <- diag(x = 10, nrow = ncol(nu0))
  }
  # Set up mhrm loop parameters
  tol <- .01
  log_lik <- NA
  #start iteration
  if (verbose_mhrm) {
    cat(
      "... burn-in completed at",
      format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
      "\n",
      sep = " "
    )
  }

  # Omega
  iter <- 0
  go <- TRUE
  if (est_omega) {

    if (verbose_mhrm) {
      cat(
        "Omega MHRM Start Time",
        format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
        "\n",
        sep = " "
      )
    }

    while (go && (iter < max_iter_mhrm)) {
      iter <- iter + 1

      # STEP 1: Stochastic imputation ------------------------------------------
      mc_draws_at_iteration_k <- mhmc_mc(
        chains = 3, y = y, obj_fun = obj_fun, link = link,
        est_omega = est_omega, est_lambda = est_lambda, est_zeta = est_zeta,
        est_nu = est_nu, omega0 = omega0, gamma0 = gamma0, lambda0 = lambda0,
        zeta0 = zeta0, nu0 = nu0, kappa0 = kappa0, omega_mu = omega_mu,
        omega_sigma2 = omega_sigma2, lambda_mu = lambda_mu,
        lambda_sigma2 = lambda_sigma2, zeta_mu = zeta_mu,
        zeta_sigma2 = zeta_sigma2, nu_mu = nu_mu, nu_sigma2 = nu_sigma2,
        burn = burn, thin = thin, min_tune = min_tune,
        tune_int = tune_int, max_tune = max_tune, niter = niter
      )

      # STEP 2: Stochastic approximation ---------------------------------------

      gain <- 1 / iter

      tmp <- deriv_omega(
        y = y,
        omega = omega0,
        omega_mu = omega_mu,
        omega_sigma2 = omega_sigma2,
        gamma = gamma0,
        lambda = if (est_lambda) {
          matrix(
            data = apply(
              X = abind::abind(
                lapply(mc_draws_at_iteration_k$mhmcdraws, function(x) {
                  x[["lambda_draws"]]
                }),
                along = 1
              ),
              MARGIN = c(2, 3),
              FUN = mean
            ),
            nrow = nrow(lambda0),
            ncol = ncol(lambda0),
            byrow = TRUE
          )
        } else {
          lambda0
        },
        zeta = if (est_zeta) {
          matrix(
            data = apply(
              X = abind::abind(
                lapply(mc_draws_at_iteration_k$mhmcdraws, function(x) {
                  x[["zeta_draws"]]
                }),
                along = 1
              ),
              MARGIN = c(2, 3),
              FUN = mean
            ),
            nrow = nrow(zeta0),
            ncol = ncol(zeta0),
            byrow = TRUE
          )
        } else {
          zeta0
        },
        zeta_mu = zeta_mu,
        zeta_sigma2 = zeta_sigma2,
        nu = matrix(data =
                      if (est_nu) {
                        array(
                          data = apply(
                            X = abind::abind(
                              lapply(
                                mc_draws_at_iteration_k$mhmcdraws, function(x) {
                                  x[["nu_draws"]]
                                }),
                              along = 1
                            ),
                            MARGIN = c(2, 3),
                            FUN = mean
                          ),
                          dim = dim(nu0)
                        )
                      } else {
                        nu0
                      },
                    ncol = 1,
                    nrow = nrow(nu0),
                    byrow = TRUE
        ),
        kappa = kappa0,
        est_zeta = est_zeta,
        link = link
      )
      grad_omega <- matrix(
        data = unlist(tmp$fpd),
        nrow = nrow(omega0),
        ncol = ncol(omega0),
        byrow = TRUE
      )
      info <- array(
        data = unlist(
          tmp$post_info),
        dim =  c(ncol(omega0), ncol(omega0), nrow(y))
      )
      info1_omega <- array(data = info0_omega, dim = dim(x = info)) +
        gain * (info - array(data = info0_omega, dim = dim(x = info)))


      # STEP 3: Robbins-Monro update -------------------------------------------
      # Note using ginv instead of solve because appears to be more stable

      # Omega
      inv_omega <- array(
        data = apply(X = info1_omega,
                     MARGIN = c(2, 3),
                     FUN = function(x) {
                       try(expr = MASS::ginv(X = x, tol = .Machine$double.xmin),
                           silent = TRUE)
                     }
        ),
        dim = dim(x = info1_omega)
      )
      if (est_omega) {
        omega1 <- NULL
        if (any(
          inv_omega == "Error in svd(X) : infinite or missing values in 'x'\n")
        ) {
          omega1 <- mhmcburn$omegaEAP
        } else {
          for (i in seq(1, nrow(x = y), 1)) {
            omega1 <- rbind(
              omega1,
              omega0[i, ] + gain * grad_omega[i, ] %*% inv_omega[, , i]
            )
          }
        }
      }


      p <- obj_fun(
        y = y,
        omega = omega1,
        gamma = gamma0,
        lambda = lambda0,
        zeta = zeta0,
        nu = nu0,
        kappa = kappa0,
        link = link
      )$p
      log_lik <- sum(x = log((p^y) * (1 - p) ^ (1 - y)), na.rm = TRUE)

      # Test for completion
      if (
        any(abs(omega1 - omega0) > tol)  &&
        iter < max_iter_mhrm
      ) {
        go <- TRUE
        #update values
        if (est_omega) {
          info0_omega <- info1_omega
          omega0 <- omega1
        }
        if (verbose_mhrm) {
          cat(
            "\r",
            "... at iteration ",
            iter,
            " logLik is ",
            format(x = round(x = log_lik, digits = 4), nsmall = 4),
            sep = ""
          )
        }
      } else if (
        all(abs(omega1 - omega0) <= tol) &&
        iter <= max_iter_mhrm
      ) {
        go <- FALSE
        if (verbose_mhrm) {
          cat(
            "\n... algorithm converged at",
            format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
            "\n",
            sep = " "
          )
        }
      } else if (iter > max_iter_mhrm) {
        if (verbose_mhrm) {
          cat(
            "\n... algorithm failed to converge at",
            format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
            "\n",
            sep = " "
          )
        }
      }
    }
  }

  # Lambda
  iter <- 0
  go <- TRUE
  if (est_lambda) {

    if (verbose_mhrm) {
      cat(
        "Lambda MHRM Start Time",
        format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
        "\n",
        sep = " "
      )
    }

    while (go && (iter < max_iter_mhrm)) {
      iter <- iter + 1

      # STEP 1: Stochastic imputation ------------------------------------------
      mc_draws_at_iteration_k <- mhmc_mc(
        chains = 3, y = y, obj_fun = obj_fun, link = link,
        est_omega = est_omega, est_lambda = est_lambda, est_zeta = est_zeta,
        est_nu = est_nu, omega0 = omega0, gamma0 = gamma0, lambda0 = lambda0,
        zeta0 = zeta0, nu0 = nu0, kappa0 = kappa0, omega_mu = omega_mu,
        omega_sigma2 = omega_sigma2, lambda_mu = lambda_mu,
        lambda_sigma2 = lambda_sigma2, zeta_mu = zeta_mu,
        zeta_sigma2 = zeta_sigma2, nu_mu = nu_mu, nu_sigma2 = nu_sigma2,
        burn = burn, thin = thin, min_tune = min_tune,
        tune_int = tune_int, max_tune = max_tune, niter = niter
      )

      # STEP 2: Stochastic approximation ---------------------------------------

      gain <- 1 / iter

      tmp <- deriv_lambda(
        y = y,
        omega = if (est_omega) {
          matrix(
            data = apply(
              X = abind::abind(
                lapply(mc_draws_at_iteration_k$mhmcdraws, function(x) {
                  x[["omega_draws"]]
                }),
                along = 1
              ),
              MARGIN = c(2, 3),
              FUN = mean
            ),
            nrow = nrow(omega0),
            ncol = ncol(omega0),
            byrow = TRUE
          )
        } else {
          omega0
        },
        gamma = gamma0,
        lambda = lambda0,
        lambda_mu = lambda_mu,
        lambda_sigma2 = lambda_sigma2,
        zeta = if (est_zeta) {
          matrix(
            data = apply(
              X = abind::abind(
                lapply(mc_draws_at_iteration_k$mhmcdraws, function(x) {
                  x[["zeta_draws"]]
                }),
                along = 1
              ),
              MARGIN = c(2, 3),
              FUN = mean
            ),
            nrow = nrow(zeta0),
            ncol = ncol(zeta0),
            byrow = TRUE
          )
        } else {
          zeta0
        },
        nu = matrix(data =
                      if (est_nu) {
                        array(
                          data = apply(
                            X = abind::abind(
                              lapply(
                                mc_draws_at_iteration_k$mhmcdraws, function(x) {
                                  x[["nu_draws"]]
                                }),
                              along = 1
                            ),
                            MARGIN = c(2, 3),
                            FUN = mean
                          ),
                          dim = dim(nu0)
                        )
                      } else {
                        nu0
                      },
                    ncol = ncol(nu0),
                    nrow = nrow(nu0),
                    byrow = TRUE
        ),
        kappa = kappa0,
        link = link
      )
      grad_lambda <- matrix(
        data = unlist(tmp$fpd),
        nrow = nrow(lambda0),
        ncol = ncol(lambda0),
        byrow = TRUE
      )
      info <- array(
        data = unlist(
          tmp$post_info),
        dim =  c(ncol(lambda0), ncol(lambda0), ncol(y))
      )
      info1_lambda <- array(data = info0_lambda, dim = dim(x = info)) +
        gain * (info - array(data = info0_lambda, dim = dim(x = info)))


      # STEP 3: Robbins-Monro update -------------------------------------------
      # Note using ginv instead of solve because appears to be more stable

      # Lambda
      inv_lambda <- array(
        data = apply(X = info1_lambda,
                     MARGIN = c(2, 3),
                     FUN = function(x) {
                       try(expr = MASS::ginv(X = x, tol = .Machine$double.xmin),
                           silent = TRUE)
                     }
        ),
        dim = dim(x = info1_lambda)
      )
      if (any(
        inv_lambda == "Error in svd(X) : infinite or missing values in 'x'\n")
      ) {
        lambda1 <- mhmcburn$lambdaEAP
      } else {
        lambda1 <- NULL
        if ((is.null(x = constraints))) {
          for (i in seq_len(length.out = ncol(x = y))) {
            lambda1 <- rbind(
              lambda1,
              lambda0[i, ] + gain * grad_lambda[i, ] %*% inv_lambda[, , i]
            )
          }
        } else {
          array_to_mat_func <- function(d) {
            nrows <- sum(sapply(X = d, FUN = NROW))
            ncols <- sum(sapply(X = d, FUN = NCOL))
            ans <- matrix(data = 0, nrow = nrows, ncol = ncols)
            i1 <- 1
            j1 <- 1
            for (m in d) {
              i2 <- i1 + NROW(m) - 1
              j2 <- j1 + NCOL(m) - 1
              ans[i1:i2, j1:j2] <- m
              i1 <- i2 + 1
              j1 <- j2 + 1
            }
            return(ans)
          }
          inv_lambda_mat <- lapply(X = seq(dim(inv_lambda)[3]),
                                   FUN = abind::asub,
                                   x = inv_lambda, dims = 3)
          inv_lambda_mat <- array_to_mat_func(d = inv_lambda_mat)
          constr_fun <- function(x) {
            out <- rep(x = NA, nrow(constraints[[3]]))
            count <- 1
            for (i in 1:nrow(constraints[[3]])) {
              for (j in 2:(sum(constraints[[3]][i, ] != 0))) {
                out[count] <- x[which(x = constraints[[3]][i, ] != 0)[1]] -
                  x[which(x = constraints[[3]][i, ] != 0)[j]]
                count <- count + 1
              }
            }
            return(out)
          }
          lambda_vec <- c(lambda0)
          A <- numDeriv::jacobian(func = constr_fun, x = lambda_vec)
          E <- lapply(seq(dim(inv_lambda)[3]), abind::asub, x = info1_lambda, dims = 3)
          E <- array_to_mat_func(E)
          E <- diag(diag(E))
          E.augment <- rbind(cbind(E, t(A)),
                             cbind(A, matrix(0, nrow(A), nrow(A))))
          E.augment.inv <- solve(E.augment)
          #E.inv <- E.augment.inv[1:ncol(A), 1:ncol(A)]
          grad_lambda <- array(data = t(c(grad_lambda)) %*% constraints[[4]],
                               dim = dim(grad_lambda))
          lambda1 <- array(
            data = c(lambda0) + gain * (E.augment.inv[1:ncol(A), 1:ncol(A)] %*%
                                          t(t(c(grad_lambda)))),
            dim = dim(lambda0)
          )
        }
      }
      p <- obj_fun(
        y = y,
        omega = omega0,
        gamma = gamma0,
        lambda = lambda1,
        zeta = zeta0,
        nu = nu0,
        kappa = kappa0,
        link = link
      )$p
      log_lik <- sum(x = log((p^y) * (1 - p) ^ (1 - y)), na.rm = TRUE)

      # Test for completion
      if (
        any(abs(lambda1 - lambda0) > tol) &&
        iter < max_iter_mhrm
      ) {
        go <- TRUE
        #update values

        if (est_lambda) {
          info0_lambda <- info1_lambda
          lambda0 <- lambda1
        }
        if (verbose_mhrm) {
          cat(
            "\r",
            "... at iteration ",
            iter,
            " logLik is ",
            format(x = round(x = log_lik, digits = 4), nsmall = 4),
            sep = ""
          )
        }
      } else if (
        all(abs(lambda1 - lambda0) <= tol) &&
        iter <= max_iter_mhrm
      ) {
        go <- FALSE
        if (verbose_mhrm) {
          cat(
            "\n... algorithm converged at",
            format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
            "\n",
            sep = " "
          )
        }
      } else if (iter > max_iter_mhrm) {
        if (verbose_mhrm) {
          cat(
            "\n... algorithm failed to converge at",
            format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
            "\n",
            sep = " "
          )
        }
      }
    }
  }

  # Nu
  iter <- 0
  go <- TRUE
  if (est_nu) {

    if (verbose_mhrm) {
      cat(
        "Nu MHRM Start Time",
        format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
        "\n",
        sep = " "
      )
    }
    while (go && (iter < max_iter_mhrm)) {
      iter <- iter + 1

      # STEP 1: Stochastic imputation ------------------------------------------
      mc_draws_at_iteration_k <- mhmc_mc(
        chains = 3, y = y, obj_fun = obj_fun, link = link,
        est_omega = est_omega, est_lambda = est_lambda, est_zeta = est_zeta,
        est_nu = est_nu, omega0 = omega0, gamma0 = gamma0, lambda0 = lambda0,
        zeta0 = zeta0, nu0 = nu0, kappa0 = kappa0, omega_mu = omega_mu,
        omega_sigma2 = omega_sigma2, lambda_mu = lambda_mu,
        lambda_sigma2 = lambda_sigma2, zeta_mu = zeta_mu,
        zeta_sigma2 = zeta_sigma2, nu_mu = nu_mu, nu_sigma2 = nu_sigma2,
        burn = burn, thin = thin, min_tune = min_tune,
        tune_int = tune_int, max_tune = max_tune, niter = niter
      )

      # STEP 2: Stochastic approximation ---------------------------------------

      gain <- 1 / iter

      tmp <- deriv_nu(
        y = y,
        omega = if (est_omega) {
          matrix(
            data = apply(
              X = abind::abind(
                lapply(mc_draws_at_iteration_k$mhmcdraws, function(x) {
                  x[["omega_draws"]]
                }),
                along = 1
              ),
              MARGIN = c(2, 3),
              FUN = mean
            ),
            nrow = nrow(omega0),
            ncol = ncol(omega0),
            byrow = TRUE
          )
        } else {
          omega0
        },
        gamma = gamma0,
        lambda = if (est_lambda) {
          matrix(
            data = apply(
              X = abind::abind(
                lapply(mc_draws_at_iteration_k$mhmcdraws, function(x) {
                  x[["lambda_draws"]]
                }),
                along = 1
              ),
              MARGIN = c(2, 3),
              FUN = mean
            ),
            nrow = nrow(lambda0),
            ncol = ncol(lambda0),
            byrow = TRUE
          )
        } else {
          lambda0
        },
        zeta = if (est_zeta) {
          matrix(
            data = apply(
              X = abind::abind(
                lapply(mc_draws_at_iteration_k$mhmcdraws, function(x) {
                  x[["zeta_draws"]]
                }),
                along = 1
              ),
              MARGIN = c(2, 3),
              FUN = mean
            ),
            nrow = nrow(zeta0),
            ncol = ncol(zeta0),
            byrow = TRUE
          )
        } else {
          zeta0
        },
        nu = nu0,
        nu_mu = nu_mu,
        nu_sigma2 = nu_sigma2,
        kappa = kappa0,
        link = link
      )
      grad_nu <- matrix(
        data = unlist(tmp$fpd),
        nrow = nrow(nu0),
        ncol = ncol(nu0),
        byrow = TRUE
      )
      info <- array(
        data = unlist(
          tmp$post_info),
        dim =  c(ncol(nu0), ncol(nu0), ncol(y))
      )
      info1_nu <- array(data = info0_nu, dim = dim(x = info)) +
        gain * (info - array(data = info0_nu, dim = dim(x = info)))


      # STEP 3: Robbins-Monro update -------------------------------------------
      # Note using ginv instead of solve because appears to be more stable

      # Nu
      inv_nu <- array(
        data = apply(X = info1_nu,
                     MARGIN = c(2, 3),
                     FUN = function(x) {
                       try(expr = MASS::ginv(X = x, tol = .Machine$double.xmin),
                           silent = TRUE)
                     }
        ),
        dim = dim(x = info1_nu)
      )
      if (est_nu) {
        nu1 <- NULL
        if (any(
          inv_nu == "Error in svd(X) : infinite or missing values in 'x'\n")
        ) {
          nu1 <- mhmcburn$nuEAP
        } else {
          if (is.null(x = constraints)) {
            for (i in 1:ncol(x = y)) {
              nu1 <- rbind(
                nu1,
                nu0[i, ] + gain * grad_nu[i, ] %*% inv_nu[, , i]
              )
            }
          } else {
            constr_fun <- function(x) {
              out <- rep(x = NA, sum(rowSums(x = constraints[[1]] != 0) - 1))
              count <- 1
              for (i in 1:nrow(constraints[[1]])) {
                if ((sum(constraints[[1]][i, ] != 0) - 1) > 0) {
                  for (j in 2:(sum(constraints[[1]][i, ] != 0))) {
                    out[count] <- x[which(x = constraints[[1]][i, ] != 0)[1], ] -  x[which(x = constraints[[1]][i, ] != 0)[j], ]
                    count <- count + 1
                  }
                }
              }
              return(out)
            }
            A <- numDeriv::jacobian(func = constr_fun, x = nu0)
            E <- array(data = 0, dim = c(nrow(x = nu0), nrow(x = nu0)))
            diag(x = E) <- info1_nu[, , 1:nrow(x = nu0)]
            E.augment <- rbind(cbind(E, t(A)),
                                cbind(A, matrix(0, nrow(A), nrow(A))))
            E.augment.inv <- solve(E.augment)
            #E.inv <- E.augment.inv[1:ncol(A), 1:ncol(A)]
            grad_nu <- array(data = t(c(grad_nu)) %*% constraints[[2]],
                             dim = dim(grad_nu))
            nu1 <- nu0 + gain * (E.augment.inv[1:ncol(A), 1:ncol(A)] %*%
                                   grad_nu)
          }
        }
      }
      p <- obj_fun(
        y = y,
        omega = omega0,
        lambda = lambda0,
        gamma = gamma0,
        zeta = zeta0,
        nu = nu1,
        kappa = kappa0,
        link = link
      )$p
      log_lik <- sum(x = log((p^y) * (1 - p) ^ (1 - y)), na.rm = TRUE)

      # Test for completion
      if (
        any(abs(nu1 - nu0) > tol) &&
        iter < max_iter_mhrm
      ) {
        go <- TRUE
        #update values
        if (est_nu) {
          info0_nu <- info1_nu
          nu0 <- nu1
        }
        if (verbose_mhrm) {
          cat(
            "\r",
            "... at iteration ",
            iter,
            " logLik is ",
            format(x = round(x = log_lik, digits = 4), nsmall = 4),
            sep = ""
          )
        }
      } else if (
        all(abs(nu1 - nu0) <= tol) &&
        iter <= max_iter_mhrm
      ) {
        go <- FALSE
        if (verbose_mhrm) {
          cat(
            "\n... algorithm converged at",
            format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
            "\n",
            sep = " "
          )
        }
      } else if (iter > max_iter_mhrm) {
        if (verbose_mhrm) {
          cat(
            "\n... algorithm failed to converge at",
            format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
            "\n",
            sep = " "
          )
        }
      }
    }
  }
  # Compute final loglikelihood value
  p <- obj_fun(
    y = y,
    omega = if (est_omega) {
      omega1
    } else {
      omega0
    },
    gamma = gamma0,
    lambda = if (est_lambda) {
      lambda1
    } else {
      lambda0
    },
    zeta = zeta0,
    nu = if (est_nu) {
      nu1
    } else {
      nu0
    },
    kappa = kappa0,
    link = link
  )$p
  log_lik <- sum(x = log(x = (p ^ y) * (1 - p) ^ (1 - y)), na.rm = TRUE)

  if (verbose_mhrm) {
    cat(
      "Final logLik is ",
      format(x = round(x = log_lik, digits = 4), nsmall = 4),
      "\n",
      sep = " "
    )
  }
  if (est_omega) {
    tmp_omega_deriv <- deriv_omega(
      y = y,
      omega = omega1,
      gamma = gamma0,
      lambda = if (est_lambda) lambda1 else lambda0,
      zeta = zeta0,
      nu = if (est_nu) nu1 else nu0,
      kappa = kappa0,
      omega_mu = omega_mu,
      omega_sigma2 = omega_sigma2,
      est_zeta = FALSE
    )
  }
  if (est_lambda) {
    tmp_lambda_deriv <- deriv_lambda(
      y = y,
      omega = if (est_omega) omega1 else omega0,
      gamma = gamma0,
      lambda = lambda1,
      zeta = zeta0,
      nu = if (est_nu) nu1 else nu0,
      kappa = kappa0,
      lambda_mu = lambda_mu,
      lambda_sigma2 = lambda_sigma2
    )
  }
  if (est_nu) {
    tmp_nu_deriv <- deriv_nu(
      y = y,
      omega = if (est_omega) omega1 else omega0,
      gamma = gamma0,
      lambda = if (est_lambda) lambda1 else lambda0,
      zeta = zeta0,
      nu = nu1,
      kappa = kappa0,
      nu_mu = nu_mu,
      nu_sigma2 = nu_sigma2
    )
  }
  return(structure(.Data = list(
    "omega1" = if (est_omega) {
      omega1
    } else {
      NULL
    },
    "info1_omega" = if (est_omega) {
      tmp_omega_deriv$post_info
    } else {
      NULL
    },
    "lambda1" = if (est_lambda) {
      lambda1
    } else {
      NULL
    },
    "info1_lambda" = if (est_lambda) {
      tmp_lambda_deriv$post_info
    } else {
      NULL
    },
    "nu1" = if (est_nu) {
      nu1
    } else {
      NULL
    },
    "info1_nu" = if (est_nu) {
      tmp_nu_deriv$post_info
    } else {
      NULL
    },
    "log_lik" = log_lik
  ),
  class = c("cog_irt"))
  )
}
