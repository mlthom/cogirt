#' Simulated Data for a Signal Detection IRT Model
#'
#' Data and parameters were simulated based on the example provided in the
#' sim_dich_response.R function.
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{y}{Matrix of dichotomous responses.}
#'   \item{ystar}{Matrix of latent response variates.}
#'   \item{nu}{Mean of the item intercept parameters (scalar).}
#'   \item{lambda}{Matrix of item structure parameters.}
#'   \item{gamma}{Matrix of experimental structure parameters.}
#'   \item{omega}{Subject-level effects of the experimental manipulation.}
#'   \item{zeta}{Condition-level prediction errors.}
#'   ...
#' }
"sdirt"
