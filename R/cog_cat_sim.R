#-------------------------------------------------------------------------------
#' Perform Simulated Computerized Adaptive Testing for Task Conditions
#'
#' This function performs simulated adapting testing using the D-optimality
#' criterion (Segall, 2009) which allows the user to focus on a subset of
#' intentional abilities (or traits).
#'
#' @param data a matrix of item responses (K by IJ). Rows should contain
#' dichotomous responses (1 or 0) for the items indexed by each column.
#' @param model an IRT model name. The options are "1p" for the one-parameter
#' model, "2p" for the two parameter model, "3p" for the three-parameter, or
#' "sdt" for a signal detection-weighted model.
#' @param guessing either a single numeric guessing value or a matrix of item
#' guessing parameters (IJ by 1). This argument is only used when model = '3p'.
#' @param contrast_codes either a matrix of experimental structure parameters
#' (JM by MN) or the name in quotes of a R stats contrast function (i.e.,
#' "contr.helmert", "contr.poly", "contr.sum", "contr.treatment", or
#' "contr.SAS"). If using the R stats contrast function items in the data matrix
#' must be arranged by condition.
#' @param num_conditions the number of conditions (required if using the R stats
#' contrast function or when constraints = TRUE).
#' @param num_contrasts the number of contrasts including intercept (required if
#' using the R stats contrast function or when constraints = TRUE).
#' @param constraints either a logical (TRUE or FALSE) indicating that item
#' parameters should be constrained to be equal over the J conditions or a 1 by
#' I vector of items that should be constrained to be equal across conditions.
#' @param key an item key vector where 1 indicates target and 2 indicates
#' distractor (IJ). Required when model = 'sdt'.
#' @param omega a matrix of true omega parameters if known. These are
#' estimated using the complete data if not supplied by the user.
#' @param item_disc a matrix of item discrimination parameters. These are
#' estimated using the complete data if not supplied by the user.
#' @param item_int a matrix of item intercept parameters. These are estimated
#' using the complete data if not supplied by the user.
#' @param conditions a list of experimental conditions that the adaptive testing
#' algorithm will choose from. The word conditions is used here to refer to a
#' single item or a group of items that are administered together prior to the
#' next iteration of adaptive testing. For cognitive experiments, multiple
#' conditions can be assigned the same experimental level (e.g., memory load
#' level).
#' @param int_par the index of the intentional parameters. That is, the column
#' of the experimental effects matrix (omega) that should be optimized.
#' @param start_conditions a vector of condition(s) that are completed prior to
#' the onset of adaptive testing.
#' @param max_conditions the maximum number of conditions to administer before
#' terminating adaptive testing. If max_conditions is specified, min_se should
#' not be.
#' @param omit_conditions a vector of conditions to be ommitted from the
#' simulation.
#' @param min_se The minimum standard error of estimate needed to terminate
#' adaptive testing. If min_see is specified, max_conditions should not be.
#' @param link the name ("logit" or "probit") of the link function to be used in
#' the model.
#'
#' @return List with elements for all parameters estimated, standard error
#' values for all parameters estimated, and the conditions selected for adaptive
#' testing.
#'
#' @references
#' Segall, D. O. (2009). Principles of Multidimensional Adaptive Testing. In W.
#' J. van der Linden & C. A. W. Glas (Eds.), \emph{Elements of Adaptive Testing}
#'  (pp. 57-75). https://doi.org/10.1007/978-0-387-85461-8_3
#'
#' @examples
#' # Adapt until minimum standard error criterion met
#' res1 <- cog_cat_sim(data = ex3$y, model = 'sdt', guessing = NULL,
#'             contrast_codes = "contr.poly", num_conditions = 10,
#'             num_contrasts = 2, constraints = NULL, key = ex3$key,
#'             omega = ex3$omega, item_disc = ex3$lambda, item_int = ex3$nu,
#'             conditions = ex3$condition, int_par = 1, start_conditions = 3,
#'             min_se = .3, link = "probit")
#'
#' # Adapt until maximum conditions criterion met
#' res2 <- cog_cat_sim(data = ex3$y, model = 'sdt', guessing = NULL,
#'             contrast_codes = "contr.poly", num_conditions = 10,
#'             num_contrasts = 2, constraints = NULL, key = ex3$key,
#'             omega = ex3$omega, item_disc = ex3$lambda, item_int = ex3$nu,
#'             conditions = ex3$condition, int_par = 1, start_conditions = 3,
#'             max_conditions = 5, link = "probit")
#'
#' # Use omit_conditions argument to restrict testing conditions (useful for
#' # setting up a non-adaptive or restricted-adaptive comparison)
#' res3 <- cog_cat_sim(data = ex3$y, model = 'sdt', guessing = NULL,
#'             contrast_codes = "contr.poly", num_conditions = 10,
#'             num_contrasts = 2, constraints = NULL, key = ex3$key,
#'             omega = ex3$omega, item_disc = ex3$lambda, item_int = ex3$nu,
#'             conditions = ex3$condition, int_par = 1, start_conditions = 3,
#'             max_conditions = 3, omit_conditions = c(2:4, 6:9),
#'             link = "probit")
#'
#' # Applied to real cpt data
#' res4 <- cog_cat_sim(data = cpt$y, model = 'sdt', guessing = NULL,
#'                    contrast_codes = "contr.poly",
#'                    num_conditions = length(unique(cpt$condition)),
#'                    num_contrasts = 2, constraints = NULL, key = cpt$key,
#'                    conditions = cpt$condition, int_par = 1,
#'                    start_conditions = 100,
#'                    max_conditions = 3, link = "probit")
#'
#' # Applied to real sopt data
#' res5 <- cog_cat_sim(data = sopt$y, model = '1p', guessing = NULL,
#'                    contrast_codes = "contr.poly",
#'                    num_conditions = length(unique(sopt$condition)),
#'                    num_contrasts = 2, constraints = NULL,
#'                    conditions = sopt$condition, int_par = 1,
#'                    start_conditions = 3,
#'                    max_conditions = 3, link = "probit")
#'
#' @export cog_cat_sim
#-------------------------------------------------------------------------------

cog_cat_sim <- function(data = NULL, model = NULL, guessing = NULL,
                        contrast_codes = NULL, num_conditions = NULL,
                        num_contrasts = NULL, constraints = NULL, key = NULL,
                        omega = NULL, item_disc = NULL, item_int = NULL,
                        conditions = NULL, int_par = NULL,
                        start_conditions = NULL, max_conditions = Inf,
                        omit_conditions = NULL, min_se = -Inf,
                        link = "probit") {

  if (is.infinite(x = max_conditions)) {
    max_conditions <- length(x = unique(x = conditions))
  }
  K <- nrow(x = data)
  if (!is.null(x = num_contrasts)) N <- num_contrasts else N <- 1
  if (model == "sdt") M <- 2 else M <- 1
  omega0 <- array(data = 0, dim = c(K, M * N))
  completed_conditions <- matrix(data = start_conditions,
                                 nrow = nrow(x = data),
                                 ncol = length(x = start_conditions),
                                 byrow = TRUE)
  iter <- 1
  crit_se <- array(data = FALSE, dim = c(K, 1))

  # STEP 1: Estimate parameters using complete data if not defined -------------

  if (is.null(x = item_int) || is.null(x = item_disc)  || is.null(x = omega)) {
    tmp_arg <- list(
      data = data, model = model, guessing = guessing,
      contrast_codes = contrast_codes, num_conditions = num_conditions,
      num_contrasts = num_contrasts, constraints = constraints, key = key,
      est_lambda = if (is.null(x = item_disc) & model %in% c("2p", "3p")) {
        TRUE
      } else {
        FALSE
      },
      est_nu = if (is.null(x = item_int)) TRUE else FALSE,
      link = link, verbose_mhrm = FALSE
    )
    tmp_arg$omega0 <- omega
    tmp_arg$lambda0 <- item_disc
    tmp_arg$nu0 <- item_int
    tmp_res <- do.call(what = cog_irt, args = tmp_arg)
    if (is.null(x = omega)) {
      omega <- tmp_res$omega1
    }
    if (is.null(x = item_int)) {
      item_int <- tmp_res$nu1
    }
    if (is.null(x = item_disc)) {
      item_disc <- tmp_res$lambda1
    }
  }

  # STEP 2: Estimate omega from starting conditions ----------------------------
  tmp_item_disc <- sweep(
    x = item_disc,
    MARGIN = 1,
    STATS = matrix(data = conditions %in% completed_conditions, ncol = 1),
    FUN = "*"
  )
  tmp_res <- cog_irt(
    data = data, model = model, guessing = guessing,
    contrast_codes = contrast_codes, num_conditions = num_conditions,
    num_contrasts = num_contrasts, constraints = constraints, key = key,
    omega0 = omega0, est_lambda = FALSE, est_nu = FALSE,
    lambda0 = tmp_item_disc, nu0 = item_int, link = link, verbose_mhrm = FALSE
  )
  se_omega <- lapply(
    X = tmp_res$info1_omega,
    FUN = function(x) {
      sqrt(x = diag(x = solve(x))[int_par])
    }
  )
  crit_se <- array(
    data = unlist(lapply(X = se_omega, FUN = function(x) all(x < min_se))),
    dim = c(K, 1)
  )
  ongoing_crit_se <- crit_se
  ongoing_omega_est <- tmp_res$omega1[, int_par, drop = FALSE]
  ongoing_se_omega <- do.call(rbind, se_omega)

  # STEP 3: Select next condition ----------------------------------------------
  while (any(!crit_se) && !(max_conditions <= iter)) {

    incomplete_conditions <- lapply(
      X = 1:K,
      FUN = function(x) {
        unique(x = conditions[!conditions %in% completed_conditions[x, ]]
        )
      }
    )
    next_condition <- array(data = NA, dim = c(K, 1))
    for (k in 1:K) {
      if (! is.null(x = omit_conditions)) {
        incomplete_conditions[[k]] <- incomplete_conditions[[k]][
          which(x = !incomplete_conditions[[k]] %in% omit_conditions)
        ]
      }
      ics <- incomplete_conditions[[k]]
      tmp_deriv <- array(data = NA, dim = c(M * N, M * N, length(ics)))
      if (!crit_se[k, ]) {
        for (ic in ics) {
          tmp_item_disc <- sweep(
            x = item_disc,
            MARGIN = 1,
            STATS = matrix(
              data = conditions %in% c(completed_conditions, ic),
              ncol = 1
            ),
            FUN = "*"
          )
          tmp_deriv[, , which(ics == ic)] <- cog_irt(
            data[k, , drop = FALSE],
            model = model,
            guessing = guessing,
            contrast_codes = contrast_codes,
            num_conditions = num_conditions,
            num_contrasts = num_contrasts,
            constraints = constraints,
            key = key,
            est_omega = FALSE,
            est_lambda = FALSE,
            est_nu = FALSE,
            est_zeta = FALSE,
            omega0 = tmp_res$omega1[k, , drop = FALSE],
            lambda0 = tmp_item_disc,
            nu0 = item_int,
            link = link,
            verbose_mhrm = FALSE
          )[["info1_omega"]][[1]]
        }
        next_condition[k, ] <- incomplete_conditions[[k]][
          which.max(x = apply(
            X = tmp_deriv[int_par, int_par, , drop = FALSE],
            MARGIN = 3,
            FUN = det
          )
          )
        ]
      }
    }

    # STEP 5: Update omega estimates and control parameters --------------------
    iter <- iter + 1
    completed_conditions <- cbind(completed_conditions, next_condition)
    tmp_item_disc <- sweep(
      x = item_disc,
      MARGIN = 1,
      STATS = matrix(data = conditions %in% completed_conditions, ncol = 1),
      FUN = "*"
    )
    tmp_res <- cog_irt(
      data = data, model = model, guessing = guessing,
      contrast_codes = contrast_codes, num_conditions = num_conditions,
      num_contrasts = num_contrasts, constraints = constraints, key = key,
      omega0 = tmp_res$omega1[, , drop = FALSE], est_lambda = FALSE,
      est_nu = FALSE, lambda0 = tmp_item_disc, nu0 = item_int, link = link,
      verbose_mhrm = FALSE
    )
    se_omega <- ifelse(
      test = crit_se,
      yes = se_omega,
      no = lapply(
        X = tmp_res$info1_omega,
        FUN = function(x)  {
          sqrt(x = diag(x = solve(x))[int_par])
        }
      )
    )
    crit_se <- array(
      data = unlist(lapply(X = se_omega, FUN = function(x) all(x < min_se))),
      dim = c(K, 1)
    )
    ongoing_omega_est <- cbind(ongoing_omega_est, tmp_res$omega1[, int_par])
    ongoing_se_omega <- cbind(ongoing_se_omega, do.call(rbind, se_omega))
    ongoing_crit_se <- cbind(ongoing_crit_se, crit_se)
  }

  return(
    structure(.Data = list(
      "model" = model,
      "omega" = omega,
      "num_conditions" = num_conditions,
      "int_par" = int_par,
      "max_conditions" = max_conditions,
      "min_se" = min_se,
      "omega1" = tmp_res$omega1,
      "info1_omega" = tmp_res$info1_omega,
      "ongoing_omega_est" = ongoing_omega_est,
      "ongoing_se_omega" = ongoing_se_omega,
      "completed_conditions" = completed_conditions
    ),
    class = c(model, "cog_cat_sim")
    )
  )

}
