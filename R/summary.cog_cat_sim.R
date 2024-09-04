#-------------------------------------------------------------------------------
#' Method of Summary for cog_cat_sim S3
#'
#' This function provides summary statistics for simulated computerized adaptive
#' testing.
#'
#' @param object An object of class 'cog_cat_sim'.
#' @param ... Additional arguments.
#'
#' @export
#-------------------------------------------------------------------------------

summary.cog_cat_sim <- function(object, ...) {
  cat(
    "\n",
    "-------------------------------------------------------------------------",
    "\n",
    "CogIRT: Results of Simulated Computerized Adaptive Testing",
    "\n",
    "-------------------------------------------------------------------------",
    "\n",
    "\n",
    "Model                            = ",
    object$model,
    "\n",
    "Number of Conditions             = ",
    object$num_conditions,
    "\n",
    "Maximum Conditions Criterion     = ",
    object$max_conditions,
    "\n",
    "Minimum Standard Error Criterion = ",
    object$min_se,
    "\n",
    "\n",
    sprintf(fmt = "%-9s", ""),
    sprintf(fmt = "%15s", "Bias"),
    sprintf(fmt = "%15s", "MAE"),
    sprintf(fmt = "%15s", "RMSE"),
    "\n",
    sep = " "
  )
  for (i in object$int_par) {

    bias <- sum(object$omega1[, i] - object$omega[, i]) /
      nrow(x = object$omega)
    mae <- sum(abs(x = object$omega1[, i] - object$omega[, i])) /
      nrow(x = object$omega)
    rmse <- sqrt(x = sum((object$omega1[, i] - object$omega[, i])^2) /
                   nrow(x = object$omega))
    cat(
      sprintf(fmt = "%-9s", paste(" Omega", i)),
      format(x = round(x = bias, digits = 3), nsmall = 3, width = 16),
      format(x = round(x = mae, digits = 3), nsmall = 3, width = 15),
      format(x = round(x = rmse, digits = 3), nsmall = 3, width = 15),
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
