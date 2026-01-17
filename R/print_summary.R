#' @title Print and Summary Methods
#' @description S3 methods for het_ensemble objects.
#' @name print_summary
#' @keywords internal
NULL

#' Print het_ensemble
#'
#' @param x A fitted \code{het_ensemble} object.
#' @param digits Number of digits to print.
#' @param ... Unused.
#' @export
print.het_ensemble <- function(x, digits = 4, ...) {
  cat("HetEffects Ensemble Model\n")
  cat("=========================\n\n")
  cat(sprintf("Outcome: %s\n", x$outcome_var))
  cat(sprintf("Treatment: %s\n", x$treatment_var))
  cat(sprintf("Covariates: %d\n", length(x$covariate_vars)))
  cat(sprintf("Learners: %d\n", length(x$learner_names)))
  cat(sprintf("Family: %s\n", x$family))
  cat(sprintf("N: %d\n\n", x$data$n))

  cat("Ensemble weights:\n")
  weights <- sort(x$weights, decreasing = TRUE)
  print(round(weights, digits = digits))
  invisible(x)
}

#' Summary het_ensemble
#'
#' @param object A fitted \code{het_ensemble} object.
#' @param digits Number of digits to print.
#' @param ... Unused.
#' @export
summary.het_ensemble <- function(object, digits = 4, ...) {
  weights <- sort(object$weights, decreasing = TRUE)
  perf <- object$cv_performance

  result <- list(
    call = object$call,
    outcome = object$outcome_var,
    treatment = object$treatment_var,
    n = object$data$n,
    family = object$family,
    learner_names = object$learner_names,
    weights = weights,
    cv_performance = perf
  )

  class(result) <- "summary.het_ensemble"
  result
}

#' Print summary.het_ensemble
#'
#' @param x A \code{summary.het_ensemble} object.
#' @param digits Number of digits to print.
#' @param ... Unused.
#' @export
print.summary.het_ensemble <- function(x, digits = 4, ...) {
  cat("HetEffects Ensemble Summary\n")
  cat("===========================\n\n")
  cat(sprintf("Outcome: %s\n", x$outcome))
  cat(sprintf("Treatment: %s\n", x$treatment))
  cat(sprintf("Family: %s\n", x$family))
  cat(sprintf("N: %d\n", x$n))
  cat(sprintf("Learners: %d\n\n", length(x$learner_names)))

  cat("Ensemble weights:\n")
  print(round(x$weights, digits = digits))
  cat("\nCross-validated performance:\n")
  print(x$cv_performance, digits = digits, row.names = FALSE)
  invisible(x)
}

#' Coefficients for het_ensemble
#'
#' @param object A fitted \code{het_ensemble} object.
#' @param ... Unused.
#' @export
coef.het_ensemble <- function(object, ...) {
  object$weights
}
