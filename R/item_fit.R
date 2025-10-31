#-------------------------------------------------------------------------------
#' Internal: compute a bundle of model-fit diagnostics
#'
#' Given a cogirt-style fitted objectect, compute several diagnostics in one call:
#' (1) item-level marginals; (2) S-X2-like conditional item fit;
#' (3) stable Q3 matrix; (4) modified Q3; (5) global SRMSR from pairwise
#' proportions; and (optionally) (6) a simple posterior predictive check.
#'
#' Everything is done inline; no separate helper functions are created.
#'
#' @param object An object of class 'cog_irt'.
#' @param n_groups Integer, target number of rest-score groups for S-X2-like check.
#' @param min_n Integer, minimum group size to include in S-X2-like check.
#' @param B_ppc Integer, number of posterior predictive replications. Use 0 to skip.
#' @param ppc_stat Character, which PPC stat to use: "tscore_var", "max_q3",
#'   or "item_means_rmsd".
#'
#' @return A list with elements:
#'   $item_marginals  data.frame
#'   $item_sx2        data.frame
#'   $q3              matrix
#'   $q3_mod          matrix
#'   $srmsr_pairs     numeric
#'   $ppc             list or NULL
#'
#' @export item_fit
#-------------------------------------------------------------------------------

item_fit <- function(
    obj,
    n_groups  = 5L,
    min_n     = 15L,
    B_ppc     = 0L,
    ppc_stat  = c("tscore_var", "max_q3", "item_means_rmsd")
) {
  ppc_stat <- match.arg(ppc_stat)

  Y <- obj$y
  N <- nrow(Y)
  I <- ncol(Y)

  # model-implied probabilities
  P <- dich_response_model(
    y      = obj$y,
    omega  = obj$omega1,
    gamma  = obj$gamma0,
    lambda = obj$lambda1,
    zeta   = obj$zeta0,
    nu     = obj$nu1,
    kappa  = obj$kappa0,
    link   = obj$link
  )$p

  ## 1) item marginals
  obs_means  <- colMeans(Y, na.rm = T)
  mod_means  <- colMeans(P, na.rm = T)
  diff_means <- obs_means - mod_means
  se_means   <- sqrt(mod_means * (1 - mod_means) / N)
  z_scores   <- diff_means / pmax(se_means, 1e-8)

  item_marginals <- data.frame(
    item = seq_len(I),
    obs  = obs_means,
    mod  = mod_means,
    diff = diff_means,
    z    = z_scores
  )

  ## 2) S-X2-like
  item_sx2_list <- vector("list", I)

  for (j in seq_len(I)) {
    # rest score for item j
    rest <- rowSums(Y, na.rm = TRUE) - Y[, j]

    # make groups by rest score
    probs <- seq(0, 1, length.out = n_groups + 1L)
    cuts  <- stats::quantile(rest, probs = probs, na.rm = TRUE)
    grp   <- cut(rest, breaks = unique(cuts), include.lowest = TRUE)
    tabs  <- split(seq_len(N), grp)

    X2 <- 0
    df <- 0

    for (g in tabs) {
      ng <- length(g)

      # skip tiny groups
      if (ng < min_n) next

      p_obs <- mean(Y[g, j], na.rm = TRUE)
      p_mod <- mean(P[g, j], na.rm = TRUE)

      var_g <- ng * p_mod * (1 - p_mod)

      # skip groups where model says prob ~0 or ~1
      if (!is.finite(var_g) || var_g < 1e-10) next

      X2 <- X2 + (ng * p_obs - ng * p_mod)^2 / var_g
      df <- df + 1L
    }

    # compute p-value only if we actually used at least 1 group
    if (df > 0) {
      pval <- stats::pchisq(X2, df = df, lower.tail = FALSE)
    } else {
      pval <- NA_real_
    }

    item_sx2_list[[j]] <- data.frame(
      item = j,
      X2   = X2,
      df   = df,
      p    = pval
    )
  }

  item_sx2 <- do.call(rbind, item_sx2_list)


  ## 3) Q3 (stable) + 4) modified Q3
  R <- (Y - P) / sqrt(pmax(P * (1 - P), 1e-6))
  R <- scale(R, center = TRUE, scale = FALSE)
  Q3 <- stats::cor(R, use = "pairwise.complete.obs")
  diag(Q3) <- NA_real_
  q3_bar  <- mean(Q3, na.rm = TRUE)
  Q3_mod  <- Q3 - q3_bar

  ## 5) global SRMSR from pairwise proportions
  obs2 <- (t(Y) %*% Y) / N
  mod2 <- (t(P) %*% P) / N
  diff2 <- obs2 - mod2
  diff2[upper.tri(diff2, diag = TRUE)] <- NA_real_
  srmsr_pairs <- sqrt(mean(diff2[!is.na(diff2)]^2))

  ## 6) posterior predictive check (auto-enable)
  ppc_obj <- NULL
  if (B_ppc <= 0L && !is.null(ppc_stat)) {
    # user changed the stat but not B -> pick a default
    B_ppc <- 200L
  }
  if (B_ppc > 0L) {
    # observed stat
    obs_stat <- switch(
      ppc_stat,
      "tscore_var" = stats::var(rowSums(Y, na.rm = T)),
      "max_q3"     = max(Q3, na.rm = TRUE),
      "item_means_rmsd" = sqrt(mean((obs_means - mod_means)^2, na.rm = TRUE))
    )

    rep_stats <- numeric(B_ppc)
    for (b in seq_len(B_ppc)) {
      Yrep <- matrix(
        stats::rbinom(N * I, size = 1L, prob = as.vector(P)),
        nrow = N, ncol = I
      )
      rep_stats[b] <- switch(
        ppc_stat,
        "tscore_var" = stats::var(rowSums(Yrep, na.rm = T)),
        "max_q3" = {
          Rr <- (Yrep - P) / sqrt(pmax(P * (1 - P), 1e-6))
          Rr <- scale(Rr, center = TRUE, scale = FALSE)
          Q3r <- stats::cor(Rr, use = "pairwise.complete.obs")
          diag(Q3r) <- NA_real_
          max(Q3r, na.rm = TRUE)
        },
        "item_means_rmsd" = {
          obs_r <- colMeans(Yrep, na.rm = T)
          sqrt(mean((obs_r - mod_means)^2, na.rm = T))
        }
      )
    }
    ppp <- mean(rep_stats >= obs_stat, na.rm = T)
    ppc_obj <- list(
      stat_name = ppc_stat,
      observed  = obs_stat,
      rep_mean  = mean(rep_stats),
      ppp       = ppp,
      rep_stats = rep_stats
    )
  }

  list(
    item_marginals = item_marginals,
    item_sx2       = item_sx2,
    q3             = Q3,
    q3_mod         = Q3_mod,
    srmsr_pairs    = srmsr_pairs,
    ppc            = ppc_obj
  )
}
