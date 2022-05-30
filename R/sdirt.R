#-------------------------------------------------------------------------------
#' Signal Detection Item Response Theory
#'
#' This function estimates the signal detection item response theory model
#' (sdirt). Users must provide a vector indicating targets versus distractors
#' (or foils). Output includes individual estimates of discriminability (dprime)
#' and centered conservative bias (ccenter). Additionally, XXXXX
#'
#-------------------------------------------------------------------------------

# sdirt <- function(y = NULL, key = NULL) {
#
#   tmp <- rmmh(
#     chains = 3,key
#     y = y],
#     obj_fun = "dich_response_model",
#     est_omega = TRUE,
#     est_nu = TRUE,
#     est_zeta = TRUE,
# #TBD    lambda = rda$lambda[which(rda$list %in% complete_lists), , drop = F],
# #TBD    kappa = rda$kappa[, which(rda$list %in% complete_lists), drop = F],
# #TBD    gamma = rda$gamma,
# #TBD    omega0 = array(data = 0, dim = c(nrow(rda$y), ncol(rda$omega_mu))),
# #TBD    nu0 = array(
# #TBD      data = 0,
# #TBD      dim = c(ncol(rda$y), 1)
# #TBD    )[which(rda$list %in% complete_lists), , drop = F],
# #TBD    zeta0 = array(data = 0, dim = c(nrow(rda$y), ncol(rda$zeta_mu))),
# #TBD    omega_mu = rda$omega_mu,
# #TBD    omega_sigma2 = rda$omega_sigma2,
# #TBD    nu_mu = matrix(rda$nu_mu),
# #TBD    nu_sigma2 = matrix(rda$nu_sigma2),
# #TBD    zeta_mu = rda$zeta_mu,
# #TBD    zeta_sigma2 = rda$zeta_sigma2,
#     burn = 0,
#     thin = 5,
#     min_tune = 0,
#     tune_int = 0,
#     max_tune = 0,
#     niter = 6,
#     verbose_rmmh = F,
#     max_iter_rmmh = 200
#   )
# }
