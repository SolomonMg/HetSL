#' Credit Claiming Experiment Data
#'
#' Experimental data from the credit claiming study in Grimmer, Messing,
#' and Westwood (2017). Includes a derived binary treatment indicator
#' \code{party_treat} for party cue vs control.
#'
#' @format A data frame with 1,074 rows and 54 columns.
#' Key variables include \code{approval}, \code{therm}, \code{cond.party},
#' \code{party_treat}, \code{pid3l}, and \code{ideo3}.
#' @source Grimmer, Messing & Westwood (2017)
"het_experiment"

#' Credit Claiming Experiment Data (Cleaned)
#'
#' Preprocessed version of \code{het_experiment} following RepCode.R, including
#' respondent filtering and derived binary outcomes used in the paper.
#'
#' @format A data frame with 1,065 rows and 55 columns.
#' Key variables include \code{approve_bi}, \code{approval}, \code{therm},
#' \code{cond.party}, \code{party_treat}, \code{pid3l}, and \code{ideo3}.
#' @source Grimmer, Messing & Westwood (2017)
"het_experiment_clean"

#' Blame Prediction Outputs from the Paper
#'
#' Precomputed predictions and ensemble weights used in the paper's
#' replication materials.
#'
#' @format A list with elements:
#' \describe{
#'   \item{Treated}{Matrix of treated predictions.}
#'   \item{Control}{Matrix of control predictions.}
#'   \item{Weights}{Named vector of ensemble weights.}
#' }
#' @source Grimmer, Messing & Westwood (2017)
"blame_preds"
