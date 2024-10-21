#-------------------------------------------------------------------------------
#' Perform Simulated Computerized Adaptive Testing
#'
#' This function performs simulated adapting testing using the D-optimality
#' criterion (Segall, 2009) which allows the user to focus on a subset of
#' intentional abilities (or traits).
#'
#' @param data A matrix of item responses (K by IJ). Rows should contain
#' dichotomous responses (1 or 0) for the items indexed by each column.
#' @param model An IRT model name. The options are "1p" for the one-parameter
#' model, "2p" for the two-parameter model, "3p" for the three-parameter model,
#' or "sdt" for a signal detection-weighted model.
#' @param guessing Either a single numeric guessing value or a matrix of item
#' guessing parameters (IJ by 1). This argument is only used when model = '3p'.
#' @param contrast_codes Either a matrix of contrast codes (JM by MN) or the
#' name in quotes of a R stats contrast function (i.e., "contr.helmert",
#' "contr.poly", "contr.sum", "contr.treatment", or "contr.SAS"). If using the R
#' stats contrast function items in the data matrix must be arranged by
#' condition.
#' @param num_conditions The total number of possible conditions (required if
#' using the R stats contrast function or when constraints = TRUE).
#' @param num_contrasts The number of contrasts, including intercept (required if
#' using the R stats contrast function or when constraints = TRUE).
#' @param constraints Either a logical (TRUE or FALSE) indicating that item
#' parameters should be constrained to be equal over the J conditions, or a 1 by
#' I vector of items that should be constrained to be equal across conditions.
#' @param key An item key vector where 1 indicates a target and 2 indicates
#' a distractor (IJ). Required when model = 'sdt'.
#' @param omega A matrix of true omega parameters if known. These are
#' estimated using the complete data if not supplied by the user.
#' @param item_disc A matrix of item discrimination parameters if known. These
#' are estimated using the complete data if not supplied by the user.
#' @param item_int A matrix of item intercept parameters if known. These are
#' estimated using the complete data if not supplied by the user.
#' @param conditions A list of experimental conditions that the adaptive testing
#' algorithm will choose from. The word "conditions" here refers to a
#' single item or a group of items that should be administered together before
#' the next iteration of adaptive testing. For cognitive experiments, multiple
#' conditions can be assigned the same experimental level (e.g., memory load
#' level).
#' @param int_par The index of the intentional parameters, i.e., the column
#' of the experimental effects matrix (omega) that should be optimized.
#' @param start_conditions A vector of condition(s) completed prior to
#' the onset of adaptive testing.
#' @param max_conditions The maximum number of conditions to administer before
#' terminating adaptive testing. If max_conditions is specified, min_se should
#' not be. Note that this is the number of additional conditions to administer
#' beyond the starting conditions.
#' @param omit_conditions A vector of conditions to be omitted from the
#' simulation.
#' @param min_se The minimum standard error of estimate needed to terminate
#' adaptive testing. If min_se is specified, max_conditions should not be.
#' @param link The name ("logit" or "probit") of the link function to be used in
#' the model.
#'
#' @return A list with elements with the model used (model), true omega
#' parameters (omega), various simulation parameters, final omega estimates
#' (omega1) and information matrices (info1_omega), ongoing estimates of omega
#' (ongoing_omega_est) and standard error of the estimates (ongoing_se_omega),
#' and completed conditions (completed_conditions).
#'
#' @references
#' Segall, D. O. (2009). Principles of Multidimensional Adaptive Testing. In W.
#' J. van der Linden & C. A. W. Glas (Eds.), \emph{Elements of Adaptive Testing}
#'  (pp. 57-75). https://doi.org/10.1007/978-0-387-85461-8_3
#'
#' @examples
#' sim_res <- cog_cat_sim(data = ex3$y, model = 'sdt', guessing = NULL,
#'                     contrast_codes = "contr.poly", num_conditions = 10,
#'                     num_contrasts = 2, constraints = NULL, key = ex3$key,
#'                     omega = ex3$omega, item_disc = ex3$lambda,
#'                     item_int = ex3$nu, conditions = ex3$condition,
#'                     int_par = c(1, 2), start_conditions = 3,
#'                     max_conditions = 3, link = "probit")
#' summary(sim_res)
#' plot(sim_res)
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
  y <- as.matrix(x = data)
  if (!is.numeric(x = y)) {
    stop("'y' contains non-numeric data.",
         call. = FALSE)
  }
  if (is.null(x = model)) {
    stop("'model' argument not specified.",
         call. = FALSE)
  }
  if (! model %in% c("1p", "2p", "3p", "sdt")) {
    stop("'model' argument is not correctly specified.",
         call. = FALSE)
  }
  if (model == "3p" && is.null(guessing)) {
    stop("'model' = '3p' but 'guessing' is NULL.",
         call. = FALSE)
  }
  if (model != "3p" && !is.null(guessing)) {
    message("'model' is not '3p'. 'guessing' argument is ignored")
    guessing <- NULL
  }
  if (!is.null(guessing) && !is.numeric(guessing)) {
    stop("'guessing' contains non-numeric data.", call. = FALSE)
  }
  if (!is.null(guessing)) {
    if ((is.matrix(guessing) && nrow(guessing) == ncol(y)) ||
        length(guessing) == 1) {
    } else {
      stop("'guessing' should either be a single numeric guessing value or a
              matrix of item guessing parameters (IJ by 1).",
           call. = FALSE)
    }
  }
  if (model %in% c("1p", "2p", "3p") && !is.null(x = key)) {
    message("'model' is not 'sdt'. 'key' argument is ignored")
  }
  if (model == "sdt") {
    if (is.null(x = key)) {
      stop("'model' sdt specified but key argument is missing.",
           call. = FALSE)
    }
    if (ncol(x = y) != length(x = key)) {
      stop("Key length does not match the number of columns in data.",
           call. = FALSE)
    }
  }
  if (model %in% c("2p", "3p")) {
    if (!is.null(x = num_conditions)) {
      if (num_conditions > 1) {
        if (is.null(x = constraints)) {
          stop("'constraints' must be TRUE when model is '2p' or '3p' and
               'num_conditions' > 1", call. = FALSE)
        } else if(is.logical(x = constraints)){
          if (constraints == FALSE) {
            stop("'constraints' must be TRUE when model is '2p' or '3p' and
               'num_conditions' > 1", call. = FALSE)
          }
        }
      }
    }
  }
  if (is.infinite(x = max_conditions)) {
    max_conditions <- length(x = unique(x = conditions))
  }
  K <- nrow(x = y)
  if (!is.null(x = num_contrasts)) N <- num_contrasts else N <- 1
  if (model == "sdt") M <- 2 else M <- 1
  omega0 <- array(data = 0, dim = c(K, M * N))
  completed_conditions <- matrix(data = start_conditions,
                                 nrow = nrow(x = y),
                                 ncol = length(x = start_conditions),
                                 byrow = TRUE)
  iter <- 1
  crit_se <- array(data = FALSE, dim = c(K, 1))

  # STEP 1: Estimate item parameters using complete data if not defined --------

  if (is.null(x = item_int) || is.null(x = item_disc)  || is.null(x = omega)) {
    cat(
      "Estimating item parameters...",
      "\n",
      sep = " "
    )
    tmp_arg <- list(
      data = y, model = model, guessing = guessing,
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
  cat(
    "Estimating omega from starting conditions...",
    "\n",
    sep = " "
  )
  tmp_item_disc <- sweep(
    x = item_disc,
    MARGIN = 1,
    STATS = matrix(data = conditions %in% completed_conditions, ncol = 1),
    FUN = "*"
  )
  tmp_res <- cog_irt(
    data = y, model = model, guessing = guessing,
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

  cat(
    "Adaptive Testing Start Time",
    format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
    "\n",
    sep = " "
  )

  # STEP 3: Select next condition ----------------------------------------------
  while (any(!crit_se) && !(max_conditions < iter)) {
    cat(
      "... CAT administration ",
      iter,
      "\n",
      sep = ""
    )
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
      data = y, model = model, guessing = guessing,
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
  cat(
    "Adaptive Testing End Time",
    format(x = Sys.time(), format = "%m/%d/%y %H:%M:%S"),
    "\n",
    sep = " "
  )

  return(
    structure(.Data = list(
      "model" = model,
      "omega" = omega,
      "num_conditions" = num_conditions,
      "int_par" = int_par,
      "start_conditions" = start_conditions,
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
