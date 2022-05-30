#-------------------------------------------------------------------------------
#' Cognitive Testing Using Computerized Adaptive Testing
#'
#' This function takes an rda file or list with select objects and returns omega
#' estimates, an information matrix, and the next best list to administer for
#' computerized adaptive testing. Adapting testing using D-optimality (see
#' Segall 2009).
#'
#' @param x .rda file (or list) containing all objects necessary to run rmmh.
#' @param obj_fun A function that calculates predictions and log-likelihood
#' values for the selected model (character).
#' @param int_par Intentional parameters. That is, the parameters to optimize
#' precision (scalar).
#'
#' @references
#' Segall, D. O. (2009). Principles of Multidimensional Adaptive Testing. In W.
#' J. van der Linden & C. A. W. Glas (Eds.), \emph{Elements of Adaptive Testing}
#'  (pp. 57-75). https://doi.org/10.1007/978-0-387-85461-8_3
#'
#' @examples
#' rda = sdirtSS
#' rda$list <- c(sapply(X = 1:(length(rda$y) / 5), FUN = rep, 5))
#' rda$y[which(!rda$list %in% c(1))] <- NA
#' cog_cat(rda = rda, obj_fun = dich_response_model, int_par = 1)
#'
#' @export cog_cat
#-------------------------------------------------------------------------------


cog_cat <- function(rda = NULL, obj_fun = NULL, int_par = NULL) {

  # STEP 1: Define complete vs. incomplete lists -------------------------------
  complete_lists <- unique(x = rda$list[which(!is.na(x = rda$y))])
  incomplete_lists <- unique(x = rda$list[which(is.na(x = rda$y))])

  # STEP 2: Estimate trait using current data ----------------------------------
  tmp <- rmmh(
    chains = 3,
    y = rda$y[1, which(rda$list %in% complete_lists), drop = F],
    obj_fun = obj_fun,
    est_omega = TRUE,
    est_nu = TRUE,
    est_zeta = TRUE,
    lambda = rda$lambda[which(rda$list %in% complete_lists), , drop = F],
    kappa = rda$kappa[, which(rda$list %in% complete_lists), drop = F],
    gamma = rda$gamma,
    omega0 = array(data = 0, dim = c(nrow(rda$y), ncol(rda$omega_mu))),
    nu0 = array(
      data = 0,
      dim = c(ncol(rda$y), 1)
    )[which(rda$list %in% complete_lists), , drop = F],
    zeta0 = array(data = 0, dim = c(nrow(rda$y), ncol(rda$zeta_mu))),
    omega_mu = rda$omega_mu,
    omega_sigma2 = rda$omega_sigma2,
    nu_mu = matrix(rda$nu_mu),
    nu_sigma2 = matrix(rda$nu_sigma2),
    zeta_mu = rda$zeta_mu,
    zeta_sigma2 = rda$zeta_sigma2,
    burn = 0,
    thin = 5,
    min_tune = 0,
    tune_int = 0,
    max_tune = 0,
    niter = 6,
    verbose_rmmh = F,
    max_iter_rmmh = 200
  )

  # STEP 3: Determine which list remaining in the bank to administer -----------
  temp <- array(
    sapply(
      incomplete_lists,
      function(x) {
        dich_response_deriv(
          y = rda$y[1, which(rda$list %in% x), drop = F],
          nu = rep(
            x = rda$nu_mu,
            length(rda$y[1, which(rda$list %in% x), drop = F])
          ),
          lambda = rda$lambda[which(rda$list %in% x), , drop = F],
          kappa = rda$kappa[, which(rda$list %in% x), drop = F],
          gamma = rda$gamma,
          omega = tmp$omega1,
          zeta = rda$zeta_mu,
          omega_mu = rda$omega_mu,
          omega_sigma2 = rda$omega_sigma2,
          zeta_mu = rda$zeta_mu,
          zeta_sigma2 = rda$zeta_sigma2
        )[["post_info"]][[1]]
      }
    ),
    dim = c(dim(rda$omega_sigma2), length(incomplete_lists))
  )

  # STEP 4: Select next list limited limit to intentional parameters -----------
  next_list <- incomplete_lists[which.max(
    x = apply(temp[int_par, int_par, , drop = F], 3, det)
  )]

  return(list(
    "omega1" = tmp$omega1,
    "info1" = tmp$info1,
    "next_list" = next_list
  ))
}
