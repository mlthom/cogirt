#-------------------------------------------------------------------------------
#' Method of Summary for cog_irt S3
#'
#' This function provides summary statistics for parameter estimates produced
#' for various CogIRT models.
#'
#' @param object An object of class 'cog_irt'.
#' @param ... Additional arguments.
#'
#' @keywords internal
#-------------------------------------------------------------------------------

summary.cog_irt <- function(object, ...) {
  if ("1p" %in% class(object)) {
    mod_name <- "One-Parameter Item Response Theory Model"
  } else if ("2p" %in% class(x = object)) {
    mod_name <- "Two-Parameter Item Response Theory Model"
  } else if ("3p" %in% class(x = object)) {
    mod_name <- "Three-Parameter Item Response Theory Model"
  } else if ("sdt" %in% class(x = object)) {
    mod_name <- "Signal Detection-Weighted IRT Model"
  } else {
    mod_name <- "User-Specified IRT Model"
  }
  cat(
    "\n",
    "-------------------------------------------------------------------------",
    "\n",
    paste("CogIRT: IRT Estimates for the ", mod_name, sep = ""),
    "\n",
    "-------------------------------------------------------------------------",
    "\n",
    "\n",
    "Number of Subjects = ",
    nrow(x = object$y),
    "\n",
    "Number of Items    = ",
    ncol(x = object$y),
    "\n",
    "log-Likelihood     = ",
    object$log_lik,
    "\n",
    "\n",
    sprintf(fmt = "%-9s", ""),
    sprintf(fmt = "%15s", "Mean"),
    sprintf(fmt = "%15s", "SD"),
    sprintf(fmt = "%15s", "Median SEE"),
    sprintf(fmt = "%15s", "Reliability"),
    "\n",
    sep = " "
  )
  for (i in 1:ncol(x = object$omega1)) {

    omega1 <- object$omega1[, i]
    mean_omega <- mean(x =  omega1)
    sd_omega <- sd(x = omega1)
    var_omega <- var(x = omega1)
    errvar_omega <- unlist(x = lapply(X = object$info1_omega,
                                      FUN = function(x)  {
                                        diag(x = solve(x))[i]
                                      }
    ))
    mean_errvar_omega <- mean(x = errvar_omega)
    se_omega <- sqrt(x = errvar_omega)
    med_se_omega <- median(x = se_omega)
    #DOI: 10.1007/S11336-011-9238-0
    rel_omega <- var_omega / (var_omega + mean_errvar_omega)
    #rel_omega <- 1 - sum(errvar_omega/var_omega) / length(x = omega1)
    cat(
      if ("sdirt" %in% class(x = object)) {
        if (i == 1) {
          sprintf(fmt = "%-9s", " D-prime")
        } else {
          sprintf(fmt = "%-9s", " C-center")
        }
      } else {
        sprintf(fmt = "%-9s", paste(" Omega", i, sep = ""))
      },
      format(x = round(x = mean_omega, digits = 3), nsmall = 3, width = 16),
      format(x = round(x = sd_omega, digits = 3), nsmall = 3, width = 15),
      format(x = round(x = med_se_omega, digits = 3), nsmall = 3, width = 15),
      format(x = round(x = rel_omega, digits = 3), nsmall = 3, width = 15),
      "\n",
      sep = " "
    )
  }
  nu <- object$nu1[, 1]
  mean_nu <- mean(x =  nu)
  sd_nu <- sd(x = nu)
  var_nu <- var(x = nu)
  errvar_nu <- unlist(x = lapply(X = object$info1_nu,
                                 FUN = function(x)  {
                                   diag(x = solve(x))[1]
                                 }
  ))
  mean_errvar_nu <- mean(x = errvar_nu)
  se_nu <- sqrt(x = errvar_nu)
  med_se_nu <- median(x = se_nu)
  (rel_nu <- var_nu / (var_nu + mean_errvar_nu))
  #(rel_nu <- 1 - sum(errvar_nu/var_nu) / length(x = nu))
  cat(
    sprintf(fmt = "%-9s", " Nu "),
    format(x = round(x = mean_nu, digits = 3), nsmall = 3, width = 16),
    format(x = round(x = sd_nu, digits = 3), nsmall = 3, width = 15),
    format(x = round(x = med_se_nu, digits = 3), nsmall = 3, width = 15),
    format(x = round(x = rel_nu, digits = 3), nsmall = 3, width = 15),
    "\n",
    sep = " "
  )
  if ("2p" %in% class(x = object) || "3p" %in% class(x = object)) {
    lambda <- object$lambda1[, 1]
    mean_lambda <- mean(x =  lambda)
    sd_lambda <- sd(x = lambda)
    var_lambda <- var(x = lambda)
    errvar_lambda <- unlist(x = lapply(X = object$info1_lambda,
                                       FUN = function(x)  {
                                         diag(x = solve(x))[1]
                                       }
    ))
    mean_errvar_lambda <- mean(x = errvar_lambda)
    se_lambda <- sqrt(x = errvar_lambda)
    med_se_lambda <- median(x = se_lambda)
    (rel_lambda <- var_lambda / (var_lambda + mean_errvar_lambda))
    #(rel_lambda <- 1 - sum(errvar_lambda/var_lambda) / length(x = lambda))
    cat(
      sprintf(fmt = "%-9s", " Lambda "),
      format(x = round(x = mean_lambda, digits = 3), nsmall = 3, width = 16),
      format(x = round(x = sd_lambda, digits = 3), nsmall = 3, width = 15),
      format(x = round(x = med_se_lambda, digits = 3), nsmall = 3, width = 15),
      format(x = round(x = rel_lambda, digits = 3), nsmall = 3, width = 15),
      "\n",
      sep = " "
    )
  }
  cat(
    "\n",
    "-------------------------------------------------------------------------",
    "\n",
    sep = " "
  )
}
