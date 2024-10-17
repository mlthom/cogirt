#-------------------------------------------------------------------------------
#' Derivatives and Information for Lambda
#'
#' This function calculates the matrix of first partial derivatives, the matrix
#' of second partial derivatives, and matrix of posterior and Fisher information
#' for the posterior distribution with respect to alpha (discrimination) based
#' on the slope-intercept form of the 1-, 2-, or 3-P item response theory model.
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
#' @param lambda_mu Mean prior for lambda (1 by JM)
#' @param lambda_sigma2 Covariance prior for lambda (JM by JM)
#' @param link Choose between "logit" or "probit" link functions.
#'
#' @return List with elements fpd (1 by JM vector of first partial derivatives
#' for alpha), spd (JM by JM matrix of second partial derivatives for alpha),
#' post_info (JM by JM posterior information matrix for alpha), and fisher_info
#' (JM by JM Fisher information matrix for alpha). Within each of these
#' elements, there are sub-elements for all IJ items
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
#'
#' @examples
#'
#'# deriv_lambda(y = ex1$y, omega = ex1$omega, gamma = ex1$gamma,
#'#              lambda = ex1$lambda, zeta = ex1$zeta, nu = ex1$nu,
#'#              lambda_mu = ex1$lambda_mu, lambda_sigma2 = ex1$lambda_sigma2,
#'#              link  = "probit")
#'
#' @section A Note About Model Notation:
#' The function converts GLLVM notation to the more typical IRT notation used by
#' Segall (1996) for ease of referencing formulas (with the exception of using
#' the slope-intercept form of the item response model).
#'
#' @keywords internal
#-------------------------------------------------------------------------------

deriv_lambda <- function(y = NULL, omega = NULL, gamma = NULL, lambda = NULL,
                         zeta = NULL, kappa = NULL, nu = NULL, lambda_mu = NULL,
                         lambda_sigma2 = NULL, link  = NULL) {
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
  # A challenge with this function is that the derivation is set up for the
  # IRT 'a' parameter. In our instance, 'a' is not a formal parameter in the
  # model but it can be defined in two ways. First, 'a' can be defined by a
  # transformation matrix that combines matrix lambda X gamma with matrix
  # lambda. This is accomplished using a transformation matrix. Then, the
  # appropriate 'theta' matrix combines matrix theta with matrix zeta. This is
  # defined as 'option1' in the code below. The problem with this option1 is
  # that the derivatives for 'a' are not equal to the derivatives for lambda.
  # But we can use the transformation matrix to covert the 'a' derivatives to
  # those for lambda, but it's not clear whether this is accurate. A second
  # option is to multiple omega into gamma and define this as 'theta'. If we do
  # this, lambda is equivalent to 'a'. This is called option 2 below. The
  # problem with this approach is that it partially ignores zeta. Zeta is used
  # to calculate p but it does not affect the scaling of the derivatives.
  # However, because zeta comprises nuisance terms, option 2 seems to be a more
  # accurate approach and used below.
  option <- 2
  if (option == 1) {
    trans_mat <- cbind(gamma, diag(1, ncol(lambda)))
    a <- lambda %*% trans_mat

    sigma2 <- diag(x = 0, nrow = ncol(lambda) + ncol(gamma))
    sigma2[seq_len(length.out = ncol(x = lambda %*% gamma)),
           seq_len(length.out = ncol(x = lambda %*% gamma))] <-
      diag(x = c(diag(lambda_sigma2) %*% abs(gamma)), nrow = ncol(gamma))
    sigma2[(1 + ncol(x = lambda %*% gamma)):ncol(x = sigma2),
           (1 + ncol(x = lambda %*% gamma)):ncol(x = sigma2)] <- lambda_sigma2
    theta <- cbind(omega, zeta)
    mu <- cbind(lambda_mu %*% gamma, lambda_mu)
  } else if (option == 2) {
    theta <- omega %*% t(gamma)
    a <- lambda
    sigma2 <- lambda_sigma2
    mu <- lambda_mu
  }
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
  for (i in seq_len(length.out = nrow(x = a))) {
    fpd[[i]] <- t(
      (
        D * apply(X = (
          theta * matrix(data = ((p[, i] - c[, i]) * (y[, i] - p[, i])),
                         nrow = nrow(x = theta),
                         ncol = ncol(x = theta),
                         byrow = FALSE)
        ) /
          matrix(data = (1 - c[, i]) * p[, i],
                 nrow = nrow(x = theta),
                 ncol = ncol(x = theta),
                 byrow = FALSE),
        MARGIN = 2,
        FUN = sum,
        na.rm = TRUE)
      ) - solve(a = sigma2) %*% t(a[i, , drop = FALSE] - mu)
    )
  }
  # Segall (1996) Equation 25; Segall (2009) Appendix; Carlson (1988) Appendix A
  spd <- list()
  for (i in seq_len(length.out = nrow(x = a))) {
    spd[[i]] <- matrix(data = NA,
                       nrow = ncol(a[i, , drop = FALSE]),
                       ncol = ncol(a[i, , drop = FALSE]))
    diag(spd[[i]]) <-  D^2 *
      apply(X =
              (
                theta^2 *
                  matrix(data = (1 - p[, i]) * ((p[, i] - c[, i]) *
                                                  (c[, i] * y[, i] - p[, i]^2)),
                         nrow = nrow(theta),
                         ncol = ncol(theta),
                         byrow = FALSE)
              ) /
              matrix(data = (p[, i]^2) * (1 - c[, i])^2,
                     nrow = nrow(theta),
                     ncol = ncol(theta),
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
  for (i in seq_len(length.out = nrow(x = a))) {
    post_info[[i]] <-
      solve(a = sigma2) +
      apply(X = D^2 *
              array(data = apply(X = theta,
                                 MARGIN = 1,
                                 FUN = function(x) {
                                   x %*% t(x)
                                 }
              ),
              dim = c(nrow(x = sigma2), nrow(x = sigma2), nrow(theta))) *
              array(data =
                      sapply(X = (1 - p[, i]) / p[, i] * ((p[, i] - c[, i]) /
                                                            (1 - c[, i]))^2,
                             FUN = rep,
                             times = nrow(x = sigma2) * nrow(x = sigma2)),
                    dim = c(nrow(x = sigma2), nrow(x = sigma2), nrow(theta))),
            MARGIN = c(1, 2),
            FUN = sum)
  }
  return(
    if (option == 1) {
      list(
        "fpd" = lapply(X = fpd, FUN = function(x) {
          x %*% MASS::ginv(trans_mat)
        }),
        "spd" = lapply(X = spd, FUN = function(x) {
          t(MASS::ginv(trans_mat)) %*% x %*% MASS::ginv(trans_mat)
        }),
        "post_info" = lapply(X = post_info, FUN = function(x) {
          t(MASS::ginv(trans_mat)) %*% x %*% MASS::ginv(trans_mat)
        }),
        "fisher_info" = lapply(X = post_info, FUN = function(x) {
          t(MASS::ginv(trans_mat)) %*% x %*% MASS::ginv(trans_mat) * -1
        })
      )
    } else if (option == 2) {
      list(
        "fpd" = fpd,
        "spd" = spd,
        "post_info" = post_info,
        "fisher_info" = lapply(X = post_info, FUN = function(x) {
          x * -1
        })
      )
    }
  )
}
