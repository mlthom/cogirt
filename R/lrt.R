#-------------------------------------------------------------------------------
#' Method of anova for cogirt S3
#'
#' This function compares fit of models produced by cogirt.
#'
#' @param object An object of class 'cogirt'.
#' @param ... Additional arguments.
#'
#' @return An object of class "anova".
#'
#' @export
#-------------------------------------------------------------------------------

lrt <- function(object, ...) {

  mcall <- match.call(expand.dots = TRUE)
  ellipsis <- list(...)
  modp <- if (length(x = ellipsis)) {
    sapply(ellipsis, inherits, "cog_irt")
  } else {
    logical(0L)
  }
  mods <- c(list(object), ellipsis[modp])
  names(mods) <- sapply(as.list(mcall)[which(c(FALSE, TRUE, modp))],
                        function(x) deparse(x))
  logLik <- sapply(mods, function(x) x$log_lik)

  logLik <- sapply(mods, function(x) x$log_lik)
  Par <- sapply(mods, function(x) x$par)
  aic <- 2 * sapply(mods, function(x) x$par) -
    2 * sapply(mods, function(x) x$log_lik)
  bic <- log(sapply(mods, function(x) nrow(x$y))) *
    sapply(mods, function(x) x$par) - 2 * sapply(mods, function(x) x$log_lik)
  chisq <- - 2 * sapply(mods, function(x) x$log_lik)
  chisq_delta <- c(NA, abs(diff(x = chisq)))
  df_delta <- c(NA, abs(diff(x = sapply(mods, function(x) x$par))))
  pvalue_delta <- c(NA, 1 - pchisq(
    q = abs(x = diff(x = chisq)),
    df = abs(x = diff(x = sapply(mods, function(x) x$par)))
  ))

  val <- data.frame(logLik = logLik, Par = Par, AIC = aic, BIC = bic,
                    "Chisq diff" = chisq_delta, "df diff" = df_delta,
                    "Pr(>Chisq)" = pvalue_delta, row.names = names(mods),
                    check.names = FALSE)

  class(val) <- c("anova", class(val))
  return(val)
}
