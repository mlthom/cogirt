#-------------------------------------------------------------------------------
#' Administer Cognitive Tests Using Computerized Adaptive Testing
#'
#' This function accepts an RDA file or a list containing selected objects and
#' returns omega estimates, the standard error of omega, and the optimal next
#' condition to administer for single-subject computerized adaptive testing.
#' Adaptive testing is guided by D-optimality (see Segall, 2009).
#'
#' @param rda An RDA file (or list) containing y, kappa, gamma, lambda,
#' condition, omega_mu, omega_sigma2, zeta_mu, zeta_sigma2, nu_mu, and
#' nu_sigma2. y should be a 1 by IJ row vector. All items not administered
#' should have NA values in y. See package documentation for definitions and
#' dimensions of these other objects.
#' @param obj_fun A function that calculates predictions and log-likelihood
#' values for the selected model (character).
#' @param int_par Intentional parameters. That is, the parameters to optimize
#' precision (scalar).
#'
#' @return A list with elements for omega parameter estimates (omega1),
#' standard error of the estimates (se_omega), and the next condition to
#' administer (next_condition).
#'
#' @references
#' Segall, D. O. (2009). Principles of Multidimensional Adaptive Testing. In W.
#' J. van der Linden & C. A. W. Glas (Eds.), \emph{Elements of Adaptive Testing}
#'  (pp. 57-75). https://doi.org/10.1007/978-0-387-85461-8_3
#'
#' @examples
#' # Define an rda file
#' rda = ex5
#' # Set all items that are not part of 3 to NA (it not already NA)
#' rda$y[which(!rda$condition %in% c(3))] <- NA
#' # Determine next condition to administer
#' cog_cat(rda = rda, obj_fun = dich_response_model, int_par = 1)
#'
#' @export cog_cat
#-------------------------------------------------------------------------------


cog_cat <- function(rda = NULL, obj_fun = NULL, int_par = NULL) {

  # STEP 1: Define complete vs. incomplete conditions --------------------------
  complete_conditions <- unique(x = rda$condition[which(!is.na(x = rda$y))])
  incomplete_conditions <- unique(x = rda$condition[which(is.na(x = rda$y))])

  # STEP 2: Estimate trait using current data ----------------------------------
  tmp_est <- mhrm(
    chains = 3,
    y = rda$y[1, which(rda$condition %in% complete_conditions), drop = FALSE],
    obj_fun = obj_fun,
    est_omega = TRUE,
    est_lambda = FALSE,
    est_zeta = TRUE,
    est_nu = TRUE,
    omega0 = array(data = 0, dim = c(nrow(rda$y), ncol(rda$omega_mu))),
    gamma0 = rda$gamma,
    lambda0 = rda$lambda[which(rda$condition %in% complete_conditions), ,
                         drop = FALSE],
    zeta0 = array(data = 0, dim = c(nrow(rda$y), ncol(rda$zeta_mu))),
    nu0 = array(
      data = 0,
      dim = c(ncol(rda$y), 1)
    )[which(rda$condition %in% complete_conditions), , drop = FALSE],
    kappa0 = rda$kappa[which(rda$condition %in% complete_conditions), ,
                       drop = FALSE],
    omega_mu = rda$omega_mu,
    omega_sigma2 = rda$omega_sigma2,
    lambda_mu = NULL,
    lambda_sigma2 = NULL,
    zeta_mu = rda$zeta_mu,
    zeta_sigma2 = rda$zeta_sigma2,
    nu_mu = rda$nu_mu,
    nu_sigma2 = rda$nu_sigma2,
    burn = 0,
    thin = 5,
    min_tune = 0,
    tune_int = 0,
    max_tune = 0,
    niter = 6,
    verbose_mhrm = FALSE,
    max_iter_mhrm = 200
  )


  # STEP 3: Determine which condition remaining in the bank to administer ------
  tmp_deriv <- array(
    sapply(
      incomplete_conditions,
      function(x) {
        deriv_omega(
          y = rda$y[1, which(rda$condition %in% x), drop = FALSE],
          omega = tmp_est$omega1,
          omega_mu = rda$omega_mu,
          omega_sigma2 = rda$omega_sigma2,
          gamma = rda$gamma,
          lambda = rda$lambda[which(rda$condition %in% x), , drop = FALSE],
          zeta = rda$zeta_mu,
          zeta_mu = rda$zeta_mu,
          zeta_sigma2 = rda$zeta_sigma2,
          nu = rep(
            x = rda$nu_mu,
            length(rda$y[1, which(rda$condition %in% x), drop = FALSE])
          ),
          kappa = rda$kappa[which(rda$condition %in% x), , drop = FALSE]
        )[["post_info"]][[1]]
      }
    ),
    dim = c(dim(rda$omega_sigma2), length(incomplete_conditions))
  )

  # STEP 4: Select next condition limited limit to intentional parameters ------
  next_condition <- incomplete_conditions[which.max(
    x = apply(tmp_deriv[int_par, int_par, , drop = FALSE], 3, det)
  )]

  return(list(
    "omega1" = tmp_est$omega1,
    "se_omega" = sqrt(x = diag(x = solve(tmp_est$info1_omega[[1]]))),
    "next_condition" = next_condition
  ))
}
