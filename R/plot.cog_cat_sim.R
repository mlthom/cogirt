#-------------------------------------------------------------------------------
#' Method of Plot for Simulated Adaptive Testing Using CogIRT S3
#'
#' This function produces plots for standard errors for
#' cog_cat_sim results
#'
#' @param x An object of class 'cog_cat_sim'.
#' @param ... Additional arguments.
#'
#' @return This function returns a base R plot displayed in the graphics device.
#' It does not return any value to the R environment.
#'
#' @export
#-------------------------------------------------------------------------------

plot.cog_cat_sim <- function(x, ...) {
  def_par <- par()["mfrow"]
  par(mfrow = n2mfrow(nr.plots = length(x = x$int_par), asp = 2))
  for (i in seq_len(length.out = length(x = x$int_par))) {
    tmp_se <- x$ongoing_se_omega[, seq(
      i,
      ncol(x = x$ongoing_se_omega),
      length(x = x$int_par)
    )]
    tmp_se <- matrix(data = unlist(x = tmp_se), ncol = ncol(x = tmp_se))
    xlim <- c(1, 1 + x$max_conditions)
    by.x <- seq(xlim[1], xlim[2], 1)
    ylim <- round(x = range(unlist(x = tmp_se)), digits = 1)
    by.y <- pretty(x = ylim, n = 7)
    plot(x = NULL, xlim = xlim, ylim = ylim, axes = FALSE,
         xlab = "Testing Stage",
         ylab = "Standard Error of Estimate",
         main = paste("omega", x$int_par[i], sep = ""))
    abline(v = by.x, col = "gray", lty = 3)
    abline(h = by.y, col = "gray", lty = 3)
    axis(side = 1, at = by.x,
         labels = c("Start", paste("CAT", 1:x$max_conditions, sep = "")),
         tick = FALSE)
    axis(side = 2, at = by.y, labels = by.y, tick = FALSE)
    for (k in seq_len(length.out = nrow(x = tmp_se))) {
      lines(x = by.x, tmp_se[k, ], col = "grey60")
    }
    lines(x = by.x, colMeans(x = tmp_se), lwd = 4)
    legend(x = xlim[2], y = ylim[2], xjust = 1,
           legend = c("Average", "Individual"), lwd = c(4, 1), bty = "n")
  }
  par(def_par)
  return(invisible(x = NULL))
}
