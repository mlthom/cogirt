#-------------------------------------------------------------------------------
#' Derivatives and Information for Nu
#'
#' This function calculates the matrix of first partial derivatives, the matrix
#' of second partial derivatives, and matrix of posterior and Fisher information
#' for the posterior distribution with respect to nu (easiness) based on the
#' slope-intercept form of the 1-, 2-, or 3-P item response theory model.
#'
#' @param y Matrix of item responses (K by IJ).
#' @param omega Examinee-level effects of the experimental manipulation (K by
#' MN).
#' @param gamma Matrix of experimental structure parameters (JM by MN).
#' @param lambda Matrix of item slope parameters (IJ by JM).
#' @param zeta Condition-level effects of the experimental manipulation (K by
#' JM).
#' @param kappa Matrix of item guessing parameters (IJ by 1). Defaults to 0.
#' @param nu Matrix of item intercept parameters (IJ by 1).
#' @param nu_mu Mean prior for nu (1 by 1)
#' @param nu_sigma2 Covariance prior for nu (1 by 1)
#' @param link Choose between "logit" or "probit" link functions.
#'
#' @return List with elements fpd (1 by 1 vector of first partial derivatives
#' for nu), spd (1 by 1 matrix of second partial derivatives for nu),
#' post_info (1 by 1 posterior information matrix for nu), and fisher_info
#' (1 by 1) Fisher information matrix for nu). Within each of these
#' elements, there are sub-elements for all IJ items.
#'
#' @section Dimensions:
#' I = Number of items per condition; J = Number of conditions; K = Number of
#' examinees; M Number of ability (or trait) dimensions; N Number of contrasts
#' (should include intercept).
#'
#' @references
#'
#' Carlson, J. E. (1987). {Multidimensional Item Response Theory Estimation: A
#' computer program} (Reprot No. ONR87-2). The American College Testing Program.
#' https://apps.dtic.mil/sti/pdfs/ADA197160.pdf
#'
#' Segall, D. O. (1996). Multidimensional adaptive testing.
#' \emph{Psychometrika, 61(2)}, 331-354. https://doi.org/10.1007/BF02294343
#'
#' Segall, D. O. (2009). Principles of Multidimensional Adaptive Testing. In W.
#' J. van der Linden & C. A. W. Glas (Eds.), \emph{Elements of Adaptive Testing}
#'  (pp. 57-75). https://doi.org/10.1007/978-0-387-85461-8_3
#'
#' @section A Note About Model Notation:
#' The function converts GLLVM notation to the more typical IRT notation used by
#' Segall (1996) for ease of referencing formulas (with the exception of using
#' the slope-intercept form of the item response model).
#'
#' @keywords internal
#-------------------------------------------------------------------------------

deriv_nu <- function(y = NULL, omega = NULL, gamma = NULL, lambda = NULL,
                     zeta = NULL, nu = NULL, kappa = NULL, nu_mu = NULL,
                     nu_sigma2 = NULL, link  = NULL) {
  link <- if (is.null(x = link)) {
    "probit"
  } else {
    link
  }
  c <- if (is.null(x = kappa)) {
    array(data = 0, dim = dim(x = y))
  } else {
    array(data = 1, dim = c(nrow(x = y), 1)) %*% t(kappa)
  }
  sigma2 <- as.matrix(nu_sigma2)
  theta <- cbind(omega, zeta)
  mu <- nu_mu
  mod <- dich_response_model(y = y, nu = nu, lambda = lambda, gamma = gamma,
                             omega = omega, zeta = zeta, kappa = kappa,
                             link  = link)
  p <- mod$p
  D <- if (link == "logit") {
    1.000
  } else if (link == "probit") {
    1.702
  }
  # Segall (1996) Equation 25; Segall (2009) Appendix; Carlson (1988) Appendix A
  fpd <- list()
  for (i in seq_len(length.out = nrow(x = nu))) {
    fpd[[i]] <- t(
      (
        D * apply(X = (
          matrix(data = ((p[, i] - c[, i]) * (y[, i] - p[, i])),
                 nrow = nrow(x = theta),
                 ncol = 1,
                 byrow = FALSE)
        ) /
          matrix(data = (1 - c[, i]) * p[, i],
                 nrow = nrow(x = theta),
                 ncol = 1,
                 byrow = FALSE),
        MARGIN = 2,
        FUN = sum,
        na.rm = TRUE)
      ) - solve(a = sigma2) %*% t(nu[i, 1, drop = FALSE] - mu)
    )
  }
  # Segall (1996) Equation 25; Segall (2009) Appendix; Carlson (1988) Appendix A
  spd <- list()
  for (i in seq_len(length.out = nrow(x = nu))) {
    spd[[i]] <- matrix(data = NA,
                       nrow = ncol(nu[i, , drop = FALSE]),
                       ncol = ncol(nu[i, , drop = FALSE]))
    diag(spd[[i]]) <-  D^2 *
      apply(X =
              (
                matrix(data = (1 - p[, i]) * ((p[, i] - c[, i]) *
                                                (c[, i] * y[, i] - p[, i]^2)),
                       nrow = nrow(x = theta),
                       ncol = 1,
                       byrow = FALSE)
              ) /
              matrix(data = (p[, i]^2) * (1 - c[, i])^2,
                     nrow = nrow(x = theta),
                     ncol = 1,
                     byrow = FALSE),
            MARGIN = 2,
            FUN = sum,
            na.rm = TRUE) - diag(solve(a = sigma2))
    spd[[i]][lower.tri(spd[[i]])] <-
      D^2 *
      apply(X = (
        apply(X = theta, MARGIN = 1, FUN = prod) *
          matrix(data = (1 - p[, i]) * ((p[, i] - c[, i]) *
                                          (c[, i] * y[, i] - p[, i]^2)),
                 nrow = nrow(theta),
                 ncol = 1,
                 byrow = FALSE)
      ) /
        matrix((p[, i]^2) * (1 - c[, i])^2,
               nrow = nrow(theta),
               ncol = 1,
               byrow = FALSE),
      MARGIN = 2,
      FUN = sum,
      na.rm = TRUE) -
      solve(a = sigma2)[lower.tri(solve(a = sigma2))]
    spd[[i]][upper.tri(spd[[i]])] <- spd[[i]][lower.tri(spd[[i]])]
  }
  # Segall (2009) Appendix Equations 3.13
  post_info <- list()
  for (i in seq_len(length.out = nrow(x = nu))) {
    post_info[[i]] <-
      solve(a = sigma2) +
      apply(X = D^2  *
              array(data = sapply(
                X = (1 - p[, i]) / p[, i] * ((p[, i] - c[, i]) /
                                               (1 - c[, i]))^2,
                FUN = rep,
                times = nrow(x = sigma2) * nrow(x = sigma2)),
                dim = c(nrow(x = sigma2), nrow(x = sigma2), nrow(theta))),
            MARGIN = c(1, 2),
            FUN = sum)
  }
  return(
    list(
      "fpd" = lapply(X = fpd, FUN = function(x) {
        x[1]
      })
      ,
      "spd" = lapply(X = spd, FUN = function(x) {
        x[1, 1]
      }),
      "post_info" = lapply(X = post_info, FUN = function(x) {
        x[1, 1]
      }),
      "fisher_info" = lapply(X = post_info, FUN = function(x) {
        x[1, 1] * -1
      })
    )
  )
}
