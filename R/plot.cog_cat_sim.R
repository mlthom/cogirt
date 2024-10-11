#-------------------------------------------------------------------------------
#' Method of Plot for Simulated Adaptive Testing Using CogIRT S3
#'
#' This function produces plots for standard errors for
#' cog_cat_sim results
#'
#' @param object An object of class 'cog_cat_sim'.
#' @param ... Additional arguments.
#'
#' @export
#-------------------------------------------------------------------------------

plot.cog_cat_sim <- function(object, ...) {
  if (object$max_conditions == 1) {
    cat("Plot is only valid for simulations with more than one iteration.")
  } else {
    def_par <- par()["mfrow"]
    par(mfrow = n2mfrow(nr.plots = length(x = object$int_par), asp = 2))
    for (i in seq_len(length.out = length(x = object$int_par))) {
      tmp_se <- object$ongoing_se_omega[, seq(
        i,
        ncol(x = object$ongoing_se_omega),
        length(x = object$int_par)
      )]
      tmp_se <- matrix(data = unlist(x = tmp_se), ncol = ncol(x = tmp_se))
      xlim <- c(1, object$max_conditions)
      by.x <- seq(xlim[1], xlim[2], 1)
      ylim <- round(x = range(unlist(x = tmp_se)), digits = 1)
      by.y <- round(x = seq(from = ylim[1],
                            to = ylim[2],
                            length.out = 7) * 10, digits = 0) / 10
      plot(x = NULL, xlim = xlim, ylim = ylim, axes = FALSE,
           xlab = "Adaptive Testing Iteration",
           ylab = "Standard Error of Estimate",
           main = paste("omega", object$int_par[i], sep = ""))
      abline(v = by.x, col = "gray", lty = 3)
      abline(h = by.y, col = "gray", lty = 3)
      axis(side = 1, at = by.x, labels = by.x, tick = FALSE)
      axis(side = 2, at = by.y, labels = by.y, tick = FALSE)
      for (k in seq_len(length.out = nrow(x = tmp_se))) {
        lines(x = by.x, tmp_se[k, ], col = "grey60")
      }
      lines(x = by.x, colMeans(x = tmp_se), lwd = 4)
      legend(x = xlim[2], y = ylim[2], xjust = 1,
             legend = c("Average", "Individual"), lwd = c(4, 1), bty = "n")
    }
    par(def_par)
  }
}
