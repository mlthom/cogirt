#-------------------------------------------------------------------------------
#' Method of Plot for cogirt S3
#'
#' This function produces plots for parameter estimates produced for various
#' cogirt models.
#'
#' @param x An x of class 'cog_irt'.
#' @param ... Additional arguments.
#'
#' @return This function returns a base R plot displayed in the graphics device.
#' It does not return any value to the R environment.
#'
#' @export
#-------------------------------------------------------------------------------

plot.cog_irt <- function(x, ...) {
  if (nrow(x = x$omega1) == 1) {
    cat("Plot is only valid for datasets with more than one subject.")
  } else {
    oldpar <- par(no.readonly = TRUE)
    on.exit(expr = par(... = oldpar))
    def_par <- par()["mfrow"]
    par(mfrow = n2mfrow(nr.plots = ncol(x = x$omega1), asp = 2))
    for (i in seq_len(ncol(x = x$omega1))) {
      omega <- x$omega1[, i]
      errvar_omega <- unlist(x = lapply(X = x$info1_omega, FUN = function(x)  {
        diag(x = solve(x))[i]
      }))
      se_omega <- sqrt(x = errvar_omega)
      xlim <- c(min(omega), max(omega))
      by.x <- pretty(x = xlim, n = 7)
      ylim <- c(0, max(se_omega))
      by.y <- pretty(x = ylim, n = 7)
      plot(x = NULL, xlim = xlim, ylim = ylim, axes = FALSE,
           xlab = paste("Omega", i, sep = ""),
           ylab = "Standard Error of Estimate", main = "")
      abline(v = by.x, col = "gray", lty = 3)
      abline(h = by.y, col = "gray", lty = 3)
      axis(side = 1, at = by.x, labels = by.x, tick = FALSE)
      axis(side = 2, at = by.y, labels = by.y, tick = FALSE)
      points(x = omega, y = se_omega, pch = 21, col = "black",
             bg = rgb(red = 0, green = 0, blue = 0, alpha = .3))
    }
    par(def_par)
  }
  return(invisible(x = NULL))
}
