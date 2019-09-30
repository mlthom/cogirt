#-------------------------------------------------------------------------------
#' MCMH Parameter Estimates for Single Chain
#'
#' This function calculates MCMH parameter estimates for a single chain. The
#' MCMH method implemented follows Patz and Junker (1999). The approach to
#' tuning the scale and covariance matrix follows BDA and SAS 9.2 User Guide,
#' 2nd Ed. "The MCMC Procedure: Tuning the Proposal Distribution".
#'
#' @references
#'
#' Patz, R. J., & Junker, B. W. (1999). A Straightforward Approach to Markov
#' Chain Monte Carlo Methods for Item Response Models. Journal of Educational
#' and Behavioral Statistics, 24(2), 146.
#'
#' @examples
#'mcmhSC(ch = NULL, y = sdirt$y, objFun = dich_response_model, estOmega = T,
#'     estNu = T, estZeta = T, lambda = sdirt$lambda, gamma = sdirt$gamma,
#'     omega0 = array(data = 0, dim = dim(sdirt$omega)),
#'     nu0 = array(data = 0, dim = c(ncol(sdirt$nu), 1)),
#'     zeta0 = array(data = 0, dim = dim(sdirt$zeta)),
#'     omega_mu = sdirt$omega_mu, omega_sigma2 = sdirt$omega_sigma2,
#'     nu_mu = matrix(sdirt$nu_mu), nu_sigma2 = matrix(sdirt$nu_sigma2),
#'     zeta_mu = sdirt$zeta_mu, zeta_sigma2 = sdirt$zeta_sigma2,
#'     burn=0, thin=10, min_tune=50, tune_int=50, max_tune=1000, niter=2000,
#'     weight=1/1, verboseMCMH=T)
#'
#' @export mcmhSC
#-------------------------------------------------------------------------------

mcmhSC <- function(ch = NULL, y = y, objFun = NULL, estOmega = T, estNu = T,
                 estZeta = T, lambda = NULL, gamma = NULL, omega0 = NULL,
                 nu0 = NULL, zeta0 = NULL, omega_mu = NULL, omega_sigma2 = NULL,
                 nu_mu = NULL, nu_sigma2 = NULL, zeta_mu = NULL,
                 zeta_sigma2 = NULL, burn = NULL, thin = NULL, min_tune = NULL,
                 tune_int = NULL, max_tune = NULL, niter = NULL, weight = NULL,
                 verboseMCMH = F) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package \"MASS\" needed for this function to work. Please install.",
         call. = FALSE)
  }
  if(verboseMCMH) {
    cat(
      "\n",
      "MCMC LOG",
      "Chain",
      ch,
      "Start Time",
      format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
      "\n\n",
      sep=" "
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
  o.scale <- 2.38 / sqrt(x = ncol(x = omega0))
  n.scale <- 2.38 / sqrt(x = ncol(x = nu0))
  z.scale <- 2.38 / sqrt(x = ncol(x = zeta0))
  # Candidate proposal variance (see BDA3 p. 296)
  o.prop <- array(
    data = diag(.75, ncol(omega0)),
    dim = c(ncol(omega0), ncol(omega0), nrow(omega0))
  )
  cand.o.var <- o.prop * o.scale
  n.prop <- array(
    data = diag(.25, 1),
    dim = c(1, 1, nrow(nu0))
  )
  cand.n.var <- n.prop * n.scale
  z.prop <- array(
    data = diag(.05, ncol(zeta0)),
    dim = c(ncol(zeta0), ncol(zeta0), nrow(zeta0))
  )
  cand.z.var <- z.prop * z.scale
  for (i in 1:niter) {
    iter_num <- i
    # Draw omegas
    if(estOmega) {
      omega1 <- omega0 + matrix(
        data = sapply(
          X = 1:nrow(omega0),
          FUN = function(i) {
            MASS::mvrnorm(
              n = 1,
              mu = rep(0, ncol(omega0)),
              Sigma = cand.o.var[, , i])
          }
        ),
        nrow = nrow(omega0),
        ncol = ncol(omega0),
        byrow = T
      )
      p0 <-objFun(y = y, nu = array(data = nu0, dim = dim(y)), lambda = lambda,
                  gamma = gamma, omega = omega0, zeta = zeta0)$p
      tmp0 <- rowSums(x = log((p0^y) * (1 - p0)^(1 - y)), na.rm = T)
      ll0 <- tmp0 + log(x = mvtnorm::dmvnorm(
        x = omega0,
        mean = omega_mu,
        sigma = omega_sigma2
      ))
      p1 <-objFun(y = y, nu = array(data = nu0, dim = dim(y)), lambda = lambda,
                  gamma = gamma, omega = omega1, zeta = zeta0)$p
      tmp1 <- rowSums(x = log((p1^y) * (1 - p1)^(1 - y)), na.rm = T)
      ll1 <- tmp1 + log(x = mvtnorm::dmvnorm(
        x = omega1,
        mean = omega_mu,
        sigma = omega_sigma2
      ))
      # Accept draws
      accept <- ll1 - ll0
      accept <-ifelse(test = accept > 0, yes = 0, no = accept)
      accept <- ifelse(
        test = runif(nrow(omega0)) < exp(x = accept),
        yes = 1,
        no = 0
      )
      omega1[which(x = accept != 1),] <- omega0[which(x = accept != 1), ]
      omegaaccept <- accept
      if(iter_num < tune_int) {
        recent_omega_draws <- abind::abind(
          recent_omega_draws,
          array(data = t(omega1), dim = c(1, ncol(omega0), nrow(omega0))),
          along = 1
        )
      } else {
        recent_omega_draws <- abind::abind(
          array(recent_omega_draws[-1,,], dim(x = recent_omega_draws) - c(1, 0, 0)),
          array(data = t(x = omega1), dim = c(1, ncol(omega0), nrow(omega0))),
          along = 1
        )
      }
      # Adjust tuning rate and candidate covariance matrix
      omegarate <- apply(
        X = recent_omega_draws,
        MARGIN = c(3),
        FUN = function(x) {mean(x = !duplicated(x = x))}
      )
      if(iter_num %in% seq(from = min_tune, to = max_tune, by = tune_int)) {
        o.scale <- (o.scale * qnorm(p = .23/2)) / qnorm(p = (omegarate * .997 + .001) /2)
        cand.o.var <- (
          ((iter_num / (max_tune + weight * max_tune))*array(
            data = apply(
              X = array(
                data = recent_omega_draws[seq(from = 1, to = tune_int,by = ceiling(x = tune_int / 10)), , ],
                dim = c(length(x = seq(from = 1, to = tune_int, by = ceiling(x = tune_int / 10))), ncol(omega0), nrow(omega0))),
              MARGIN = c(3),
              FUN = var),
            dim = c(ncol(omega0), ncol(omega0), nrow(omega0)))) +
            ((1 - (iter_num / (max_tune + weight * max_tune)))) * array(
              data = unlist(lapply(X = o.scale, function(x) {diag(x = x, ncol(omega0))})),
              dim = c(ncol(omega0), ncol(omega0), nrow(omega0)))
        )
      }
      if(iter_num > burn & iter_num %in% seq(from = burn + 1, to = niter, by = thin)) {
        omega_draws <- abind::abind(
          omega_draws,
          array(data = t(omega1), dim = c(1, ncol(omega0), nrow(omega0))),
          along = 1)
      }
    } else {
      omega1 <- omega0
    }
    # Draw nus
    if(estNu) {
      nu1 <- nu0 + matrix(
        data = sapply(
          X = 1:nrow(nu0),
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
      p0 <-objFun(y = y, nu = array(data = nu0, dim = dim(y)), lambda = lambda,
                  gamma = gamma, omega = omega0, zeta = zeta0)$p
      tmp0 <- colSums(x = log((p0^y) * (1 - p0)^(1 - y)), na.rm = T)
      ll0 <- tmp0 + log(x = mvtnorm::dmvnorm(
        x = nu0,
        mean = nu_mu,
        sigma = nu_sigma2
      ))
      p1 <-objFun(y = y, nu = array(data = nu1, dim = dim(y)), lambda = lambda,
                  gamma = gamma, omega = omega0, zeta = zeta0)$p
      tmp1 <- colSums(x = log((p1 ^ y) * (1 - p1)^(1 - y)), na.rm = T)
      ll1 <- tmp1 + log(x = mvtnorm::dmvnorm(
        x = nu1,
        mean = nu_mu,
        sigma = nu_sigma2
      ))
      # Accept draws
      accept <- ll1 - ll0
      accept <-ifelse(test = accept > 0, yes = 0, no = accept)
      accept <- ifelse(test = runif(nrow(nu0)) < exp(x = accept), yes = 1, no = 0)
      nu1[which(x = accept != 1), ] <- nu0[which(x = accept != 1), ]
      nuaccept <- accept
      if(iter_num < tune_int) {
        recent_nu_draws <- abind::abind(
          recent_nu_draws,
          array(data = t(nu1), dim = c(1, 1, nrow(nu0))),
          along = 1
        )
      } else {
        recent_nu_draws <- abind::abind(
          array(data = recent_nu_draws[-1, , ], dim = dim(x = recent_nu_draws) - c(1, 0, 0)),
          array(data = t(nu1), dim = c(1, 1, nrow(nu0))),
          along=1
        )
      }
      # Adjust tuning rate and candidate covariance matrix
      nurate <- apply(
        X = recent_nu_draws,
        MARGIN = c(3),
        FUN = function(x) {mean(x = !duplicated(x))}
      )
      if(iter_num %in% seq(from = min_tune, to = max_tune, by = tune_int)) {
        n.scale <- (n.scale * qnorm(p = .23 / 2)) / qnorm(p = (nurate * .997 + .001) / 2)
        cand.n.var <- (
          ((iter_num / (max_tune + weight * max_tune))*array(
            data = apply(
              X = array(data = recent_nu_draws[seq(from = 1, to = tune_int,by = ceiling(x = tune_int / 10)), , ], dim = c(length(x = seq(from = 1, to = tune_int, by = ceiling(x = tune_int / 10))), ncol(nu0), nrow(nu0))),
              MARGIN = c(3),
              FUN = var),
            dim = c(ncol(nu0), ncol(nu0), nrow(nu0))
          )) +
            ((1 - (iter_num / (max_tune + weight * max_tune)))) * array(
              data = unlist(lapply(X = n.scale, function(x) {diag(x = x, ncol(nu0))})),
              dim = c(ncol(nu0), ncol(nu0), nrow(nu0)))
        )
      }
      if(iter_num > burn & iter_num %in% seq(from = burn+1, to = niter, by = thin)) {
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
    if(estZeta) {
      zeta1 <- zeta0 + matrix(
        data = sapply(
          X = 1:nrow(zeta0),
          FUN = function(i) {MASS::mvrnorm(
            n = 1,
            mu = rep(0, ncol(zeta0)),
            Sigma = cand.z.var[, , i])
          }
        ),
        nrow = nrow(zeta0),
        ncol = ncol(zeta0),
        byrow = T
      )
      p0 <-objFun(y = y, nu = array(data = nu0, dim = dim(y)), lambda = lambda,
                  gamma = gamma, omega = omega0, zeta = zeta0)$p
      tmp0 <- rowSums(x = log((p0^y) * (1 - p0)^(1 - y)), na.rm = T)
      ll0 <- tmp0 + log(x = mvtnorm::dmvnorm(
        x = zeta0,
        mean = zeta_mu,
        sigma = zeta_sigma2
      ))
      p1 <-objFun(y = y, nu = array(data = nu0, dim = dim(y)), lambda = lambda,
                  gamma = gamma, omega = omega0, zeta = zeta1)$p
      tmp1 <- rowSums(x = log((p1^y) * (1 - p1)^(1 - y)), na.rm = T)
      ll1 <- tmp1 + log(x = mvtnorm::dmvnorm(
        x = zeta1,
        mean = zeta_mu,
        sigma = zeta_sigma2
      ))
      # Accept draws
      accept <- ll1 - ll0
      accept <-ifelse(test = accept > 0, yes = 0, no = accept)
      accept <- ifelse(
        test = runif(nrow(zeta0)) < exp(x = accept),
        yes = 1,
        no = 0
      )
      zeta1[which(x = accept != 1),] <- zeta0[which(x = accept != 1), ]
      zetaaccept <- accept
      if(iter_num < tune_int) {
        recent_zeta_draws <- abind::abind(
          recent_zeta_draws,
          array(data = t(zeta1), c(1, ncol(zeta0), nrow(zeta0))),
          along = 1
        )
      } else {
        recent_zeta_draws <- abind::abind(
          array(recent_zeta_draws[-1,,], dim(x = recent_zeta_draws)-c(1,0,0)),
          array(t(x = zeta1), c(1, ncol(zeta0), nrow(zeta0))),
          along = 1
        )
      }
      # Adjust tuning rate and candidate covariance matrix
      zetarate <- apply(
        X = recent_zeta_draws,
        MARGIN = c(3),
        FUN = function(x) {mean(x = !duplicated(x = x))}
      )
      if(iter_num %in% seq(from = min_tune, to = max_tune, by = tune_int)) {
        z.scale <- (z.scale * qnorm(p = .23/2)) / qnorm(p = (zetarate * .997 + .001) /2)
        cand.z.var <- (
          ((iter_num / (max_tune + weight * max_tune))*array(
            data = apply(X = array(data = recent_zeta_draws[seq(from = 1, to = tune_int,by = ceiling(x = tune_int / 10)), , ],dim = c(length(x = seq(from = 1, to = tune_int, by = ceiling(x = tune_int / 10))), ncol(zeta0), nrow(zeta0))),
                         MARGIN = c(3),
                         FUN = var),
            dim = c(ncol(zeta0), ncol(zeta0), nrow(zeta0))
          )) +
            ((1 - (iter_num / (max_tune + weight * max_tune)))) * array(
              data = unlist(lapply(X = z.scale, FUN = function(x) {diag(x = x, ncol(zeta0))})),
              dim = c(ncol(zeta0), ncol(zeta0), nrow(zeta0)))
        )
      }
      if(iter_num > burn & iter_num %in% seq(from = burn + 1, to = niter, by = thin)) {
        zeta_draws <- abind::abind(
          zeta_draws, array(data = t(zeta1), dim = c(1, ncol(zeta0), nrow(zeta0))),
          along = 1
        )
      }
    } else {
      zeta1 <- zeta0
    }
    #Update logLikelihood and estimates
    ll <- objFun(y = y, nu = array(data = nu1, dim = dim(y)), lambda = lambda,
                 gamma = gamma, omega = omega1, zeta = zeta1)$loglikelihood
    omega0 <- omega1
    nu0 <- nu1
    zeta0 <- zeta1
    if(i %in% seq(from = 1, to = niter, by = thin) & verboseMCMH) {
      draws <- round(x = c(omega0, nu0, zeta0, ll), digits = 4)
      cat("Iter.",
          format(x = i, width = 4, justify = "left"), sep="",
          append = TRUE
      )
      cat(" Ch.",
          format(x = ch, width = 2),
          sep="",
          append = TRUE
      )
      cat(" ll =",
          format(x = round(x = ll, digits = 2), nsmall = 2, width = 10),
          sep = "",
          append=TRUE
      )
      if(estOmega) {
        cat(" Mean omega Rate = ",
            format(x = round(x = 100 * mean(x = omegarate)), width = 3),
            "% ",
            "Mean o.scale = ",
            format(x = round(x = mean(x = o.scale), 3),width = 6, nsmall = 3),
            sep="",
            append=TRUE
        )
      }
      if(estNu) {
        cat(" Mean nu Rate = ",
            format(x = round(x = 100 * mean(x = nurate)), width = 3),
            "% ",
            "Mean n.scale = ",
            format(x = round(x = mean(x = n.scale), digits = 3),
                   width = 6,
                   nsmall = 3
            ),
            sep="",
            append=TRUE
        )
      }
      if(estZeta) {
        cat(" Mean zeta Rate = ",
            format(x = round(x = 100 * mean(x = zetarate)), width = 3),
            "% ",
            "Mean z.scale = ",
            format(x = round(x = mean(x = z.scale), digits = 3),
                   width = 6,
                   nsmall =3
            ),
            sep="",
            append=TRUE)
      }
      cat(" Time = ",
          format(x = Sys.time(), format = "%H:%M:%S"),
          sep="",
          append=TRUE
      )
      cat("\n", sep=",", append=TRUE)
    }
  }
  return(list(
    "omega_draws" = omega_draws, "nu_draws" = nu_draws,
    "zeta_draws" = zeta_draws, "cand.o.var" = cand.o.var,
    "cand.n.var" = cand.n.var, "cand.z.var" = cand.z.var
  ))
}
