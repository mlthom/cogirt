#' Simulated Data for a Unidimensional Two-Parameter Item Response Model
#'
#' Data and parameters were simulated based on example 1 provided for the
#' sim_dich_response.R function.
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{y}{Matrix of dichotomous responses.}
#'   \item{ystar}{Matrix of latent response variates.}
#'   \item{omega}{Subject-level effects of the experimental manipulation.}
#'   \item{omega_mu}{Vector of means for the subject-level effects of the
#'   experimental manipulation (1 by K * M).}
#'   \item{omega_sigma2}{Covariance matrix for the subject-level effects of the
#'   experimental manipulation (K * M by K * M).}
#'   \item{gamma}{Matrix of experimental structure parameters.}
#'   \item{lambda}{Matrix of item structure parameters.}
#'   \item{lambda_mu}{Vector of means for the item structure parameters
#'   (1 by JM).}
#'   \item{lambda_sigma2}{Covariance matrix for the item structure parameters
#'   (JM by JM).}
#'   \item{nu}{Mean of the item intercept parameters (scalar).}
#'   \item{nu_mu}{Mean of the item intercept parameters (scalar).}
#'   \item{nu_sigma2}{Variance of the item intercept parameters (scalar).}
#'   \item{zeta}{Condition-level prediction errors.}
#'   \item{zeta_mu}{Vector of means for the condition-level prediction errors
#'   (1 by J * M).}
#'   \item{zeta_sigma2}{Covariance matrix for the condition-level prediction
#'   errors (J * M by J * M).}
#'   \item{kappa}{Matrix of item guessing parameters (K by IJ).}
#'   \item{condition}{Condition vector indiciting distinct conditions or time
#'   points.}
#'   \item{key}{Item key vector where 1 indicates target and 2 indicates
#'   distractor (IJ)}
#'   ...
#' }
"ex1"

#' Simulated Data for a Signal Detection Weighted IRT Model
#'
#' Data and parameters were simulated based on example 2 provided for the
#' sim_dich_response.R function.
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{y}{Matrix of dichotomous responses.}
#'   \item{ystar}{Matrix of latent response variates.}
#'   \item{omega}{Subject-level effects of the experimental manipulation.}
#'   \item{omega_mu}{Vector of means for the subject-level effects of the
#'   experimental manipulation (1 by K * M).}
#'   \item{omega_sigma2}{Covariance matrix for the subject-level effects of the
#'   experimental manipulation (K * M by K * M).}
#'   \item{gamma}{Matrix of experimental structure parameters.}
#'   \item{lambda}{Matrix of item structure parameters.}
#'   \item{lambda_mu}{Vector of means for the item structure parameters
#'   (1 by JM).}
#'   \item{lambda_sigma2}{Covariance matrix for the item structure parameters
#'   (JM by JM).}
#'   \item{nu}{Mean of the item intercept parameters (scalar).}
#'   \item{nu_mu}{Mean of the item intercept parameters (scalar).}
#'   \item{nu_sigma2}{Variance of the item intercept parameters (scalar).}
#'   \item{zeta}{Condition-level prediction errors.}
#'   \item{zeta_mu}{Vector of means for the condition-level prediction errors
#'   (1 by J * M).}
#'   \item{zeta_sigma2}{Covariance matrix for the condition-level prediction
#'   errors (J * M by J * M).}
#'   \item{kappa}{Matrix of item guessing parameters (K by IJ).}
#'   \item{condition}{Condition vector indiciting distinct conditions or time
#'   points.}
#'   \item{key}{Item key vector where 1 indicates target and 2 indicates
#'   distractor (IJ)}
#'   ...
#' }
"ex2"

#' Simulated Data for a Signal Detection Weighted IRT Model with an Experimental
#' Design
#'
#' Data and parameters were simulated based on example 3 provided for the
#' sim_dich_response.R function.
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{y}{Matrix of dichotomous responses.}
#'   \item{ystar}{Matrix of latent response variates.}
#'   \item{omega}{Subject-level effects of the experimental manipulation.}
#'   \item{omega_mu}{Vector of means for the subject-level effects of the
#'   experimental manipulation (1 by K * M).}
#'   \item{omega_sigma2}{Covariance matrix for the subject-level effects of the
#'   experimental manipulation (K * M by K * M).}
#'   \item{gamma}{Matrix of experimental structure parameters.}
#'   \item{lambda}{Matrix of item structure parameters.}
#'   \item{lambda_mu}{Vector of means for the item structure parameters
#'   (1 by JM).}
#'   \item{lambda_sigma2}{Covariance matrix for the item structure parameters
#'   (JM by JM).}
#'   \item{nu}{Mean of the item intercept parameters (scalar).}
#'   \item{nu_mu}{Mean of the item intercept parameters (scalar).}
#'   \item{nu_sigma2}{Variance of the item intercept parameters (scalar).}
#'   \item{zeta}{Condition-level prediction errors.}
#'   \item{zeta_mu}{Vector of means for the condition-level prediction errors
#'   (1 by J * M).}
#'   \item{zeta_sigma2}{Covariance matrix for the condition-level prediction
#'   errors (J * M by J * M).}
#'   \item{kappa}{Matrix of item guessing parameters (K by IJ).}
#'   \item{condition}{Condition vector indiciting distinct conditions or time
#'   points.}
#'   \item{key}{Item key vector where 1 indicates target and 2 indicates
#'   distractor (IJ)}
#'   ...
#' }
"ex3"

#' Simulated Data for a Unidimensional Two-Parameter Item Response Model
#' with Two Measurement Occasions
#'
#' Data and parameters were simulated based on example 4 provided for the
#' sim_dich_response.R function.
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{y}{Matrix of dichotomous responses.}
#'   \item{ystar}{Matrix of latent response variates.}
#'   \item{omega}{Subject-level effects of the experimental manipulation.}
#'   \item{omega_mu}{Vector of means for the subject-level effects of the
#'   experimental manipulation (1 by K * M).}
#'   \item{omega_sigma2}{Covariance matrix for the subject-level effects of the
#'   experimental manipulation (K * M by K * M).}
#'   \item{gamma}{Matrix of experimental structure parameters.}
#'   \item{lambda}{Matrix of item structure parameters.}
#'   \item{lambda_mu}{Vector of means for the item structure parameters
#'   (1 by JM).}
#'   \item{lambda_sigma2}{Covariance matrix for the item structure parameters
#'   (JM by JM).}
#'   \item{nu}{Mean of the item intercept parameters (scalar).}
#'   \item{nu_mu}{Mean of the item intercept parameters (scalar).}
#'   \item{nu_sigma2}{Variance of the item intercept parameters (scalar).}
#'   \item{zeta}{Condition-level prediction errors.}
#'   \item{zeta_mu}{Vector of means for the condition-level prediction errors
#'   (1 by J * M).}
#'   \item{zeta_sigma2}{Covariance matrix for the condition-level prediction
#'   errors (J * M by J * M).}
#'   \item{kappa}{Matrix of item guessing parameters (K by IJ).}
#'   \item{condition}{Condition vector indiciting distinct conditions or time
#'   points.}
#'   \item{key}{Item key vector where 1 indicates target and 2 indicates
#'   distractor (IJ)}
#'   ...
#' }
"ex4"


#' Simulated Single Subject Data for a Signal Detection Weighted IRT Model with
#' an Experimental Design
#'
#' Data and parameters were simulated based on example 5 provided for the
#' sim_dich_response.R function.
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{y}{Matrix of dichotomous responses.}
#'   \item{ystar}{Matrix of latent response variates.}
#'   \item{omega}{Subject-level effects of the experimental manipulation.}
#'   \item{omega_mu}{Vector of means for the subject-level effects of the
#'   experimental manipulation (1 by K * M).}
#'   \item{omega_sigma2}{Covariance matrix for the subject-level effects of the
#'   experimental manipulation (K * M by K * M).}
#'   \item{gamma}{Matrix of experimental structure parameters.}
#'   \item{lambda}{Matrix of item structure parameters.}
#'   \item{lambda_mu}{Vector of means for the item structure parameters
#'   (1 by JM).}
#'   \item{lambda_sigma2}{Covariance matrix for the item structure parameters
#'   (JM by JM).}
#'   \item{nu}{Mean of the item intercept parameters (scalar).}
#'   \item{nu_mu}{Mean of the item intercept parameters (scalar).}
#'   \item{nu_sigma2}{Variance of the item intercept parameters (scalar).}
#'   \item{zeta}{Condition-level prediction errors.}
#'   \item{zeta_mu}{Vector of means for the condition-level prediction errors
#'   (1 by J * M).}
#'   \item{zeta_sigma2}{Covariance matrix for the condition-level prediction
#'   errors (J * M by J * M).}
#'   \item{kappa}{Matrix of item guessing parameters (K by IJ).}
#'   \item{condition}{Condition vector indiciting distinct conditions or time
#'   points.}
#'   \item{key}{Item key vector where 1 indicates target and 2 indicates
#'   distractor (IJ)}
#'   ...
#' }
"ex5"

#' N-Back Data
#'
#' N-Back task accuracy data collected from an online experiment. The condition
#' vector indicates working memory load level (1-back, 2-back, 3-back, or
#' 4-back). The key indicates whether items are targets (1) or distractors (2).
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{y}{Matrix of dichotomous responses.}
#'   \item{key}{Item key vector where 1 indicates target and 2 indicates
#'   distractor (IJ)}
#'   \item{condition}{Condition vector indiciting distinct conditions or time
#'   points.}
#' }
"nback"

#' Sternberg Data
#'
#' Sternberg task accuracy data collected from an online experiment. The
#' condition vector indicates working memory load level (2, 4, 6, 8, 10, or 12
#' items).The key indicates whether items are targets (1) or distractors (2).
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{y}{Matrix of dichotomous responses.}
#'   \item{key}{Item key vector where 1 indicates target and 2 indicates
#'   distractor (IJ)}
#'   \item{condition}{Condition vector indiciting distinct conditions or time
#'   points.}
#' }
"sternberg"

#' CPT Data
#'
#' CPT task accuracy data collected from an online experiment. The
#' condition vector indicates backward mask onset (50, 100, 150, or 200 ms).The
#' key indicates whether items are targets (1) or distractors (2).
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{y}{Matrix of dichotomous responses.}
#'   \item{key}{Item key vector where 1 indicates target and 2 indicates
#'   distractor (IJ)}
#'   \item{condition}{Condition vector indiciting distinct conditions or time
#'   points.}
#' }
"cpt"

#' Flanker Data
#'
#' Flanker task accuracy data collected from an online experiment. The
#' condition vector indicates level of congruency ("congruent, incongruent_part,
#' incongruent_all, neutral).
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{y}{Matrix of dichotomous responses.}
#'   \item{condition}{Condition vector indiciting distinct conditions or time
#'   points.}
#' }
"flanker"

#' SOPT Data
#'
#' Self-Ordered Pointing Task (SOPT) accuracy data collected from an online
#' experiment. The condition vector indicates working memory load level (3, 6,
#' 9, or 12 items).
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{y}{Matrix of dichotomous responses.}
#'   \item{condition}{Condition vector indiciting distinct conditions or time
#'   points.}
#' }
"sopt"

#' PLT Data
#'
#' Probabilistic Learning Task (SOPT) accuracy data collected from an online
#' experiment. The condition vector indicates feedback consistency (90%, 80%,
#' 70%, 60%). The targ vector indicates which side is the target item. The fdbk
#' vector indicates which side was rewarded.
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{y}{Matrix of dichotomous responses.}
#'   \item{targ}{Item targ left vs. right vector (IJ)}
#'   \item{fdbk}{Item fdbk left vs. right vector (IJ)}
#'   \item{condition}{Condition vector indiciting distinct conditions or time
#'   points.}
#' }
"plt"
