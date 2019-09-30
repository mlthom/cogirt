#-------------------------------------------------------------------------------
#' Derivatives and Information for the Dichotomous Response Model
#'
#' This function calculates the matrix of first partial dervitives, matrix of
#' second partial dervitives, and information matrix for the posterior
#' distribution with respect to omega The formulas are based on Segall (1996;
#' 2009).
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
#' dich_response_deriv(y = sdirt$y, nu = sdirt$nu, lambda = sdirt$lambda,
#'                 gamma = sdirt$gamma, omega = sdirt$omega, zeta = sdirt$zeta,
#'                 omega_mu = sdirt$omega_mu, omega_sigma2 = sdirt$omega_sigma2,
#'                 zeta_mu = sdirt$zeta_mu, zeta_sigma2 = sdirt$zeta_sigma2)
#'
#' @section A Note About Model Notation:
#' The function converts GLLVM notation to the more typical IRT notation used by
#' Segall (1996) for ease of referencing formulas.

#' @export dich_response_deriv
#-------------------------------------------------------------------------------

dich_response_deriv <- function(y, nu, lambda, gamma, omega, zeta, omega_mu,
                                omega_sigma2, zeta_mu, zeta_sigma2,
                                link  = 'logit') {
  a <- cbind(lambda %*% gamma, lambda)
  ro <- nrow(omega_sigma2)
  rz <- nrow(zeta_sigma2)
  ra <- nrow(a)
  ca <- ncol(a)
  sigma2 <- diag(x = 0, nrow = ro + rz)
  sigma2[1:ro, 1:ro] <- omega_sigma2
  sigma2[(1 + ro):(ro + rz), (1 + ro):(ro + rz)] <- zeta_sigma2
  rs <- nrow(sigma2)
  theta <- cbind(omega, zeta)
  mu <- c(omega_mu, zeta_mu)
  mod <- dich_response_model(y = y, nu = nu, lambda = lambda, gamma = gamma,
                             omega = omega, zeta = zeta, link  = 'logit')
  p <- mod$p
  # Note that this is the OPPOSITE of the simulation. I.e., if the variance of
  # the latent response variate was 1.702 (logit) in the simualtion, then no
  # adjustment is needed for a logit model (i.e., D = 1.000). Conversely, if the
  # varaince was 1.000 (probit), than a 1.702 adjustment is needed.
  D <- if(link == 'logit') {
    1.000
  } else if(link == 'probit'){
    1.702
  }
  # Segall (1996) Equation 25; Segall (2009) Appendix
  fpd <- list()
  for(i in 1:nrow(theta)){
    fpd[[i]] <- t(
      (
        D * apply(X = (
          a * matrix(data = (p[i, ] * (y[i, ] - p[i, ])),
                     nrow = ra,
                     ncol = ca,
                     byrow = F)
        ) /
          matrix(data = p[i, ],
                 nrow = ra,
                 ncol = ca,
                 byrow = F),
        MARGIN = 2,
        FUN = sum)
      ) - solve(sigma2) %*% t(theta[i, , drop = F] - mu)
    )
  }
  # Segall (1996) Equations 30 & 31; Segall (2009) Appendix
  spd <- list()
  for(i in 1:nrow(theta)){
    spd[[i]] <- matrix(data = NA,
                       nrow = ncol(theta[i, , drop = F]),
                       ncol = ncol(theta[i, , drop = F]))
    diag(spd[[i]]) <-  D^2 *
      apply(X =
              (
                a^2 *
                  matrix(data = (1 - p[i, ]) * (p[i, ] * (-p[i, ]^2)),
                         nrow = nrow(a),
                         ncol = ncol(a),
                         byrow = F)
              )/
              matrix(data = p[i, ]^2,
                     nrow = nrow(a),
                     ncol = ncol(a),
                     byrow = F),
            MARGIN = 2,
            FUN = sum) - diag(solve(sigma2))
    spd[[i]][lower.tri(spd[[i]])] <-
      D^2 *
      apply(X = (
        apply(X = a, MARGIN = 1, FUN = prod) *
          matrix(data = (1-p[i, ]) * (p[i, ] * (-p[i, ]^2)),
                 nrow = nrow(a),
                 ncol = 1,
                 byrow = F)
      ) /
        matrix(p[i, ]^2,
               nrow = nrow(a),
               ncol = 1,
               byrow = F),
      MARGIN = 2,
      FUN = sum) -
      solve(sigma2)[lower.tri(solve(sigma2))]
    spd[[i]][upper.tri(spd[[i]])] <- spd[[i]][lower.tri(spd[[i]])]
  }
  # Segall (2009) Appendix Equations 3.13
  post_info <- list()
  for(i in 1:nrow(theta)){
    post_info[[i]] <-
      solve(sigma2) +
      apply(X = D^2 *
              array(data = apply(X = a,
                                 MARGIN = 1,
                                 FUN = function(x){x %*% t(x)}),
                    dim = c(rs, rs, nrow(a))) *
              array(data = sapply(X = (1 - p) / p * p^2,
                                  FUN = rep,
                                  times = rs*rs),
                    dim = c(rs, rs, nrow(a))),
            MARGIN = c(1, 2),
            FUN = sum)
  }
  # Return estimates for lambda only
  return(list(fpd = lapply(X = fpd, FUN = function(x){x[1:ro]}),
              spd = lapply(X = spd, FUN = function(x){x[1:ro, 1:ro]}),
              post_info = lapply(X = post_info, FUN = function(x){x[1:ro, 1:ro]}),
              fisher_info = lapply(X = post_info, FUN = function(x){x[1:ro, 1:ro] * -1})
  ))
}

# # Should be deleted before realease. Used for testing.
# test_dervis_function <- function(y, nu, lambda, gamma, omega, zeta, omega_mu,
#                                  omega_sigma2, zeta_mu, zeta_sigma2,
#                                  link  = 'logit'){
#   # First create a sum likelihood function that depends only on theta (combined
#   # omega and zeta). Then use my function then compare with numDeriv.
#   a <- cbind(lambda %*% gamma, lambda)
#   ro <- nrow(omega_sigma2)
#   rz <- nrow(zeta_sigma2)
#   ra <- nrow(a)
#   ca <- ncol(a)
#   sigma2 <- diag(x = 0, nrow = ro + rz)
#   sigma2[1:ro, 1:ro] <- omega_sigma2
#   sigma2[(1 + ro):(ro + rz), (1 + ro):(ro + rz)] <- zeta_sigma2
#   rs <- nrow(sigma2)
#   theta <- cbind(omega, zeta)
#   mu <- c(omega_mu, zeta_mu)
#   mod <- dich_response_model(y = y, nu = nu, lambda = lambda, gamma = gamma,
#                              omega = omega, zeta = zeta, link  = 'logit')
#   p <- mod$p
#   tmp_ll <- function(theta, y, nu, a) {
#     yhatstar <- nu + theta %*% t(a)
#     p <- if(link == 'logit') {
#       plogis(yhatstar)
#     } else if(link == 'probit'){
#       pnorm(yhatstar)
#     }
#     ll <- sum(log((p^y) * (1 - p)^(1 - y)))
#     return(ll)
#   }
#   fpd <- list()
#   tmp_ll(theta[1, ], y[1, ], nu[1, ], a)
#   for(i in 1:nrow(theta)){
#     fpd[[i]] <- numDeriv::jacobian(func = tmp_ll, theta[i, ], y=y[i, ],
#                                    nu=nu[i, ], a=a)
#   }
#   tmp_ll(theta[1, ], y[1, ], nu[1, ], a)
#   spd <- list()
#   for(i in 1:nrow(theta)){
#     spd[[i]] <- numDeriv::hessian(func = tmp_ll, theta[i, ], y=y[i, ],
#                                   nu=nu[i, ], a=a)
#   }
#   return(list(fpd = lapply(X = fpd, FUN = function(x){x[1:ro]}),
#     spd = lapply(X = spd, FUN = function(x){x[1:ro, 1:ro]})))
# }
#
# test_derivs <- function(index){
#   tmp1 <- dich_response_deriv(y = sdirt$y, nu = sdirt$nu, lambda = sdirt$lambda,
#                               gamma = sdirt$gamma, omega = sdirt$omega,
#                               zeta = sdirt$zeta, omega_mu = sdirt$omega_mu,
#                               omega_sigma2 = diag(10000, 6),
#                               zeta_mu = sdirt$zeta_mu,
#                               zeta_sigma2 = diag(10000, 10))
#
#   tmp2 <- test_dervis_function(y = sdirt$y, nu = sdirt$nu, lambda = sdirt$lambda,
#                                gamma = sdirt$gamma, omega = sdirt$omega,
#                                zeta = sdirt$zeta, omega_mu = sdirt$omega_mu,
#                                omega_sigma2 = diag(10000, 6),
#                                zeta_mu = sdirt$zeta_mu,
#                                zeta_sigma2 = diag(10000, 10))
#
#   par(mfrow = c(1, 2))
#   p.plot <- plot(
#     x = unlist(lapply(1:nrow(sdirt$y),function(x){tmp1[["fpd"]][[x]][index]})),
#     y = unlist(lapply(1:nrow(sdirt$y),function(x){tmp2[["fpd"]][[x]][index]})),
#     main = "fpd",
#     xlab = "CogIRT",
#     ylab = "numDeriv"
#   )
#   abline(a=0, b=1)
#   p.plot <- plot(
#     x = unlist(lapply(1:nrow(sdirt$y),function(x){
#       tmp1[["spd"]][[x]][index, index]
#     })),
#     y = unlist(lapply(1:nrow(sdirt$y),function(x){
#       tmp2[["spd"]][[x]][index, index]
#     })),
#     main = "spd",
#     xlab = "CogIRT",
#     ylab = "numDeriv"
#   )
#   abline(a=0, b=1)
#   return(p.plot)
# }
