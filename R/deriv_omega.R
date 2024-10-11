#-------------------------------------------------------------------------------
#' Derivatives and Information for Omega
#'
#' This function calculates the matrix of first partial derivatives, the matrix
#' of second partial derivatives, and matrix of posterior and Fisher information
#' for the posterior distribution with respect to omega (ability) based on the
#' slope-intercept form of the 1-, 2-, or 3-parameter item response theory
#' model.
#'
#' @param y Matrix of item responses (K by IJ).
#' @param omega Examinee-level effects of the experimental manipulation (K by
#' MN).
#' @param gamma Matrix of experimental structure parameters (JM by MN).
#' @param lambda Matrix of item structure parameters (IJ by JM).
#' @param zeta Condition-level effects of the experimental manipulation (K by
#' JM).
#' @param nu Matrix of item intercept parameters (IJ by 1).
#' @param kappa Matrix of item guessing parameters (IJ by 1). Defaults to 0.
#' @param omega_mu Mean prior for omega (1 by MN).
#' @param omega_sigma2 Covariance prior for omega (MN by MN).
#' @param zeta_mu Mean prior for zeta (1 by JM).
#' @param zeta_sigma2 Covariance prior for zeta (JM by JM).
#' @param est_zeta Logical indicating whether or not to estimate zeta
#' derivatives
#' @param link Choose between "logit" or "probit" link functions.
#'
#' @return List with elements fpd (1 by MN vector of first partial derivatives
#' for omega), spd (MN by MN matrix of second partial derivatives for omega),
#' post_info (MN by MN posterior information matrix for omega), and fisher_info
#' (MN by MN Fisher information matrix for omega). Within each of these
#' elements, there are sub-elements for all K examinees.
#'
#' @section Dimensions:
#' I = Number of items per condition; J = Number of conditions; K = Number of
#' examinees; M Number of ability (or trait) dimensions; N Number of contrasts
#' (should include intercept).
#'
#' @references
#'
#' Segall, D. O. (1996). Multidimensional adaptive testing.
#' \emph{Psychometrika, 61(2)}, 331-354. https://doi.org/10.1007/BF02294343
#'
#' Segall, D. O. (2009). Principles of Multidimensional Adaptive Testing. In W.
#' J. van der Linden & C. A. W. Glas (Eds.), \emph{Elements of Adaptive Testing}
#'  (pp. 57-75). https://doi.org/10.1007/978-0-387-85461-8_3
#'
#' @examples
#'# deriv_omega(y = ex1$y, omega = ex1$omega, gamma = ex1$gamma,
#'#             lambda = ex1$lambda, zeta = ex1$zeta, nu = ex1$nu,
#'#             omega_mu = ex1$omega_mu, omega_sigma2 = ex1$omega_sigma2,
#'#             zeta_mu = ex1$zeta_mu, zeta_sigma2 = ex1$zeta_sigma2,
#'#             est_zeta = FALSE, link  = "probit")
#'
#' @section A Note About Model Notation:
#' The function converts GLLVM notation to the more typical IRT notation used by
#' Segall (1996) for ease of referencing formulas (with the exception of using
#' the slope-intercept form of the item response model).
#'
#' @keywords internal
#-------------------------------------------------------------------------------

deriv_omega <- function(y = NULL, omega = NULL, gamma = NULL, lambda = NULL,
                        zeta = NULL, nu = NULL, kappa = NULL, omega_mu = NULL,
                        omega_sigma2 = NULL, zeta_mu = NULL, zeta_sigma2 = NULL,
                        est_zeta = TRUE, link  = NULL) {
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

  if (est_zeta) {
    a <- cbind(lambda %*% gamma, lambda)
    sigma2 <- diag(x = 0, nrow = nrow(x = omega_sigma2) + nrow(zeta_sigma2))
    sigma2[seq_len(nrow(x = omega_sigma2)), seq_len(nrow(x = omega_sigma2))] <-
      omega_sigma2
    sigma2[(1 + nrow(x = omega_sigma2)):(nrow(x = omega_sigma2)
                                         + nrow(zeta_sigma2)),
           (1 + nrow(x = omega_sigma2)):(nrow(x = omega_sigma2)
                                         + nrow(zeta_sigma2))] <- zeta_sigma2
    theta <- cbind(omega, zeta)
    mu <- cbind(omega_mu, zeta_mu)
  } else {
    a <- lambda %*% gamma
    sigma2 <- diag(x = 0, nrow = nrow(x = omega_sigma2))
    sigma2[seq_len(nrow(x = omega_sigma2)), seq_len(nrow(x = omega_sigma2))] <-
      omega_sigma2
    theta <- omega
    mu <- omega_mu
  }

  mod <- dich_response_model(y = y, nu = nu, lambda = lambda, gamma = gamma,
                             omega = omega, zeta = zeta, link  = link)
  p <- mod$p
  D <- if (link == "logit") {
    1.000
  } else if (link == "probit") {
    1.702
  }
  # Segall (1996) Equation 25; Segall (2009) Appendix
  fpd <- list()
  for (i in seq_len(length.out = nrow(x = theta))) {
    fpd[[i]] <- t(
      (
        D * apply(X = (
          a * matrix(data = ((p[i, ] - c[i, ]) * (y[i, ] - p[i, ])),
                     nrow = nrow(x = a),
                     ncol = ncol(x = a),
                     byrow = FALSE)
        ) /
          matrix(data = (1 - c[i, ]) * p[i, ],
                 nrow = nrow(x = a),
                 ncol = ncol(x = a),
                 byrow = FALSE),
        MARGIN = 2,
        FUN = sum,
        na.rm = TRUE)
      ) - solve(a = sigma2) %*% t(theta[i, , drop = FALSE] - mu)
    )
  }
  # Segall (1996) Equations 30 & 31; Segall (2009) Appendix
  spd <- list()
  for (i in seq_len(length.out = nrow(x = theta))) {
    spd[[i]] <- matrix(data = NA,
                       nrow = ncol(x = theta[i, , drop = FALSE]),
                       ncol = ncol(x = theta[i, , drop = FALSE]))
    diag(spd[[i]]) <-  D^2 *
      apply(X =
              (
                a^2 *
                  matrix(data = (1 - p[i, ]) * ((p[i, ] - c[i, ]) *
                                                  (c[i, ] * y[i, ] - p[i, ]^2)),
                         nrow = nrow(x = a),
                         ncol = ncol(x = a),
                         byrow = FALSE)
              ) /
              matrix(data = (p[i, ]^2) * (1 - c[i, ])^2,
                     nrow = nrow(x = a),
                     ncol = ncol(x = a),
                     byrow = FALSE),
            MARGIN = 2,
            FUN = sum,
            na.rm = TRUE) - diag(x = solve(a = sigma2))
    spd[[i]][lower.tri(spd[[i]])] <-
      D^2 *
      apply(X = (
        apply(X = a, MARGIN = 1, FUN = prod) *
          matrix(data = (1 - p[i, ]) * ((p[i, ] - c[i, ]) *
                                          (c[i, ] * y[i, ] - p[i, ]^2)),
                 nrow = nrow(x = a),
                 ncol = 1,
                 byrow = FALSE)
      ) /
        matrix((p[i, ]^2) * (1 - c[i, ])^2,
               nrow = nrow(x = a),
               ncol = 1,
               byrow = FALSE),
      MARGIN = 2,
      FUN = sum,
      na.rm = TRUE) -
      solve(a = sigma2)[lower.tri(solve(a = sigma2))]
    spd[[i]][upper.tri(x = spd[[i]])] <- spd[[i]][lower.tri(x = spd[[i]])]
  }
  # Segall (2009) Appendix Equations 3.13
  post_info <- list()
  for (i in seq_len(length.out = nrow(x = theta))) {
    post_info[[i]] <-
      solve(a = sigma2) +
      apply(X = D^2 *
              array(data = apply(X = a,
                                 MARGIN = 1,
                                 FUN = function(x) {
                                   x %*% t(x)
                                 }
              ),
              dim = c(nrow(x = sigma2), nrow(x = sigma2), nrow(x = a))) *
              array(data = sapply(X = (1 - p[i, ]) / p[i, ] *
                                    ((p[i, ] - c[i, ]) / (1 - c[i, ]))^2,
                                  FUN = rep,
                                  times = nrow(x = sigma2) * nrow(x = sigma2)),
                    dim = c(nrow(x = sigma2), nrow(x = sigma2), nrow(x = a))),
            MARGIN = c(1, 2),
            FUN = sum)
  }
  return(
    list(
      "fpd" = lapply(X = fpd, FUN = function(x) {
        x[seq_len(nrow(x = omega_sigma2))]
      }),
      "spd" = lapply(X = spd, FUN = function(x) {
        x[seq_len(nrow(x = omega_sigma2)), seq_len(nrow(x = omega_sigma2))]
      }),
      "post_info" = lapply(X = post_info, FUN = function(x) {
        x[seq_len(nrow(x = omega_sigma2)), seq_len(nrow(x = omega_sigma2))]
      }),
      "fisher_info" = lapply(X = post_info, FUN = function(x) {
        x[seq_len(nrow(x = omega_sigma2)), seq_len(nrow(x = omega_sigma2))] * -1
      })
    )
  )
}
