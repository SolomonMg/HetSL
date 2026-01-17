#' @title Treatment Effect Estimation
#' @description Functions for computing average and heterogeneous treatment effects.
#' @name treatment_effects
NULL

#' Average Treatment Effect
#'
#' Computes the average treatment effect (ATE) from a fitted ensemble.
#'
#' @param object A fitted \code{het_ensemble} object.
#' @param newdata Optional data frame for prediction. If NULL, uses training data.
#' @param control Value indicating the control condition (default 0).
#' @param treated Value indicating the treated condition (default 1).
#' @param se Logical. Compute standard error (default TRUE).
#' @param ci_level Confidence level for interval (default 0.95).
#' @param weights Optional observation weights for ATE calculation.
#'
#' @return A list of class \code{ate} with components:
#' \describe{
#'   \item{estimate}{Point estimate of ATE}
#'   \item{se}{Standard error (if se = TRUE)}
#'   \item{ci_lower}{Lower confidence bound}
#'   \item{ci_upper}{Upper confidence bound}
#'   \item{ci_level}{Confidence level}
#'   \item{n}{Sample size}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- het_ensemble(y ~ x1 + x2 | treat, data = mydata)
#' ate(fit)
#' }
ate <- function(object, newdata = NULL, control = 0, treated = 1,
                se = TRUE, ci_level = 0.95, weights = NULL) {

  if (!inherits(object, "het_ensemble")) {
    stop("'object' must be a het_ensemble object.", call. = FALSE)
  }
  if (is.matrix(object$data$treatment)) {
    stop("ATE is only supported for binary treatments.", call. = FALSE)
  }

  # Get CATE for all observations
  cate_values <- cate(object, newdata = newdata, control = control,
                      treated = treated)

  # Use provided weights or defaults
  if (is.null(weights)) {
    if (is.null(newdata)) {
      weights <- object$data$weights
    }
    if (is.null(weights)) {
      weights <- rep(1, length(cate_values))
    }
  }

  # Compute weighted average
  estimate <- .weighted_mean(cate_values, weights)
  n <- length(cate_values)

  result <- list(
    estimate = estimate,
    n = n,
    cate = cate_values,
    weights = weights
  )

  if (se) {
    # Analytic SE via influence function
    se_result <- .ate_se(object, cate_values, weights)
    result$se <- se_result$se
    result$influence <- se_result$influence

    # Confidence interval
    z <- qnorm(1 - (1 - ci_level) / 2)
    result$ci_lower <- estimate - z * result$se
    result$ci_upper <- estimate + z * result$se
    result$ci_level <- ci_level
  }

  class(result) <- "ate"
  result
}

#' Conditional Average Treatment Effect
#'
#' Computes the conditional average treatment effect (CATE) for each observation.
#'
#' @param object A fitted \code{het_ensemble} object.
#' @param newdata Optional data frame for prediction. If NULL, uses training data.
#' @param control Value indicating the control condition (default 0).
#' @param treated Value indicating the treated condition (default 1).
#'
#' @return Numeric vector of CATE estimates, one per observation.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- het_ensemble(y ~ x1 + x2 | treat, data = mydata)
#' cate_est <- cate(fit)
#' }
cate <- function(object, newdata = NULL, control = 0, treated = 1) {
  if (!inherits(object, "het_ensemble")) {
    stop("'object' must be a het_ensemble object.", call. = FALSE)
  }
  if (is.matrix(object$data$treatment)) {
    stop("CATE is only supported for binary treatments.", call. = FALSE)
  }

  if (is.null(newdata)) {
    X <- object$data$X
  } else {
    X <- newdata[, object$covariate_vars, drop = FALSE]
  }

  # Get counterfactual predictions
  cf <- .predict_counterfactual(object, X)
  cf$cate
}

#' Marginal Conditional Average Treatment Effect
#'
#' Computes the marginal conditional average treatment effect (MCATE)
#' across values of a moderator variable.
#'
#' @param object A fitted \code{het_ensemble} object.
#' @param moderator Character. Name of the moderator variable.
#' @param moderator_values Optional vector of moderator values at which to
#'   compute MCATE. If NULL, uses observed values (continuous) or all levels (factor).
#' @param newdata Optional data frame. If NULL, uses training data.
#' @param control Value indicating the control condition (default 0).
#' @param treated Value indicating the treated condition (default 1).
#' @param se Logical. Compute standard errors (default TRUE).
#' @param ci_level Confidence level (default 0.95).
#'
#' @return A data frame of class \code{mcate} with columns:
#' \describe{
#'   \item{moderator_value}{Values of the moderator}
#'   \item{estimate}{MCATE point estimates}
#'   \item{se}{Standard errors (if se = TRUE)}
#'   \item{ci_lower}{Lower confidence bounds}
#'   \item{ci_upper}{Upper confidence bounds}
#'   \item{n}{Number of observations at each value}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- het_ensemble(y ~ age + female + education | treat, data = mydata)
#' mcate(fit, moderator = "age")
#' mcate(fit, moderator = "female")
#' }
mcate <- function(object, moderator, moderator_values = NULL, newdata = NULL,
                  control = 0, treated = 1, se = TRUE, ci_level = 0.95) {

  if (!inherits(object, "het_ensemble")) {
    stop("'object' must be a het_ensemble object.", call. = FALSE)
  }
  if (is.matrix(object$data$treatment)) {
    stop("MCATE is only supported for binary treatments.", call. = FALSE)
  }

  # Get data
  if (is.null(newdata)) {
    data <- object$data$X
    weights <- object$data$weights
  } else {
    data <- newdata
    weights <- NULL
  }

  # Check moderator exists
  if (!(moderator %in% names(data))) {
    stop(sprintf("Moderator '%s' not found in data.", moderator), call. = FALSE)
  }

  mod_var <- data[[moderator]]

  # Get CATE values
  cate_values <- cate(object, newdata = newdata, control = control, treated = treated)

  # Determine moderator values
  if (is.null(moderator_values)) {
    if (is.factor(mod_var) || is.character(mod_var)) {
      moderator_values <- unique(mod_var)
    } else {
      # Use quantiles for continuous
      moderator_values <- quantile(mod_var, probs = seq(0, 1, by = 0.1), na.rm = TRUE)
      moderator_values <- unique(moderator_values)
    }
  }

  # Compute MCATE at each moderator value
  results <- lapply(moderator_values, function(mv) {
    if (is.factor(mod_var) || is.character(mod_var)) {
      idx <- which(mod_var == mv)
    } else {
      # For continuous, use local averaging (closest values)
      # Or bin-based averaging
      idx <- which(abs(mod_var - mv) <= diff(range(mod_var, na.rm = TRUE)) * 0.05)
      if (length(idx) < 5) {
        # Fall back to k nearest
        diffs <- abs(mod_var - mv)
        idx <- order(diffs)[1:min(20, length(diffs))]
      }
    }

    if (length(idx) == 0) {
      return(data.frame(
        moderator_value = mv,
        estimate = NA_real_,
        se = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        n = 0L
      ))
    }

    w <- if (!is.null(weights)) weights[idx] else rep(1, length(idx))
    est <- .weighted_mean(cate_values[idx], w)
    n_obs <- length(idx)

    if (se && n_obs > 1) {
      # SE via weighted variance
      var_est <- .weighted_var(cate_values[idx], w) / n_obs
      se_est <- sqrt(var_est)
      z <- qnorm(1 - (1 - ci_level) / 2)
      ci_lo <- est - z * se_est
      ci_hi <- est + z * se_est
    } else {
      se_est <- NA_real_
      ci_lo <- NA_real_
      ci_hi <- NA_real_
    }

    data.frame(
      moderator_value = mv,
      estimate = est,
      se = se_est,
      ci_lower = ci_lo,
      ci_upper = ci_hi,
      n = n_obs
    )
  })

  result <- do.call(rbind, results)
  result$moderator <- moderator

  class(result) <- c("mcate", "data.frame")
  attr(result, "moderator") <- moderator
  attr(result, "ci_level") <- ci_level

  result
}

#' Marginal Average Treatment Effect by Treatment Factor
#'
#' For multi-valued treatments, computes the average effect of each
#' treatment level relative to control.
#'
#' @param object A fitted \code{het_ensemble} object.
#' @param treatment_factor Character. Name of the multi-valued treatment variable.
#' @param control Value of the control condition.
#' @param se Logical. Compute standard errors (default TRUE).
#' @param ci_level Confidence level (default 0.95).
#'
#' @return A data frame of class \code{mate} with treatment-level effects.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # For multi-arm trials
#' mate(fit, treatment_factor = "arm")
#' }
mate <- function(object, treatment_factor = NULL, control = NULL,
                 se = TRUE, ci_level = 0.95) {

  if (!inherits(object, "het_ensemble")) {
    stop("'object' must be a het_ensemble object.", call. = FALSE)
  }

  # Get treatment from object if not specified
  if (is.null(treatment_factor)) {
    treat <- object$data$treatment_original
  } else {
    treat <- object$data$X[[treatment_factor]]
  }

  if (is.null(treat)) {
    stop("Cannot find treatment variable.", call. = FALSE)
  }

  # Get unique treatment levels
  levels <- unique(treat)

  if (is.null(control)) {
    # Use most common value or smallest as control
    if (is.numeric(levels)) {
      control <- min(levels)
    } else {
      control <- names(sort(table(treat), decreasing = TRUE))[1]
    }
  }

  treated_levels <- setdiff(levels, control)

  # Get CATE values
  cate_values <- cate(object)
  weights <- object$data$weights
  if (is.null(weights)) weights <- rep(1, length(cate_values))

  # Compute effect for each treatment level
  results <- lapply(treated_levels, function(tl) {
    # Find observations that received this treatment or control
    # MATE is average effect among those who received treatment
    idx <- which(treat == tl)

    if (length(idx) == 0) {
      return(data.frame(
        treatment = tl,
        control = control,
        estimate = NA_real_,
        se = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        n = 0L
      ))
    }

    w <- weights[idx]
    est <- .weighted_mean(cate_values[idx], w)
    n_obs <- length(idx)

    if (se && n_obs > 1) {
      var_est <- .weighted_var(cate_values[idx], w) / n_obs
      se_est <- sqrt(var_est)
      z <- qnorm(1 - (1 - ci_level) / 2)
      ci_lo <- est - z * se_est
      ci_hi <- est + z * se_est
    } else {
      se_est <- NA_real_
      ci_lo <- NA_real_
      ci_hi <- NA_real_
    }

    data.frame(
      treatment = tl,
      control = control,
      estimate = est,
      se = se_est,
      ci_lower = ci_lo,
      ci_upper = ci_hi,
      n = n_obs
    )
  })

  result <- do.call(rbind, results)
  class(result) <- c("mate", "data.frame")
  attr(result, "ci_level") <- ci_level

  result
}

#' Convenience Wrapper for Treatment Effects
#'
#' Computes multiple treatment effect summaries in one call.
#'
#' @param object A fitted \code{het_ensemble} object.
#' @param moderators Optional character vector of moderator names for MCATE.
#' @param se Logical. Compute standard errors.
#' @param ci_level Confidence level.
#'
#' @return A list of class \code{treatment_effects} containing:
#' \describe{
#'   \item{ate}{Average treatment effect}
#'   \item{cate}{Conditional average treatment effects}
#'   \item{mcate}{Marginal conditional effects (if moderators specified)}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' effects <- treatment_effects(fit, moderators = c("age", "female"))
#' }
treatment_effects <- function(object, moderators = NULL, se = TRUE,
                              ci_level = 0.95) {

  if (!inherits(object, "het_ensemble")) {
    stop("'object' must be a het_ensemble object.", call. = FALSE)
  }

  result <- list(
    ate = ate(object, se = se, ci_level = ci_level),
    cate = cate(object)
  )

  if (!is.null(moderators)) {
    result$mcate <- lapply(moderators, function(m) {
      mcate(object, moderator = m, se = se, ci_level = ci_level)
    })
    names(result$mcate) <- moderators
  }

  class(result) <- "treatment_effects"
  result
}

#' Print method for ate
#'
#' @param x A fitted \code{ate} object.
#' @param digits Number of digits to print.
#' @param ... Unused.
#' @export
print.ate <- function(x, digits = 4, ...) {
  cat("Average Treatment Effect (ATE)\n")
  cat("==============================\n\n")
  cat(sprintf("Estimate: %.*f\n", digits, x$estimate))

  if (!is.null(x$se)) {
    cat(sprintf("Std. Error: %.*f\n", digits, x$se))
    cat(sprintf("%d%% CI: [%.*f, %.*f]\n",
                round(x$ci_level * 100), digits, x$ci_lower, digits, x$ci_upper))
  }

  cat(sprintf("N: %d\n", x$n))
  invisible(x)
}

#' Print method for mcate
#'
#' @param x A fitted \code{mcate} object.
#' @param digits Number of digits to print.
#' @param ... Unused.
#' @export
print.mcate <- function(x, digits = 4, ...) {
  mod <- attr(x, "moderator")
  cat(sprintf("Marginal Conditional ATE by %s\n", mod))
  cat(paste(rep("=", 40), collapse = ""), "\n\n")

  # Format for printing without recursive print dispatch
  x_print <- x[, c("moderator_value", "estimate", "se", "ci_lower", "ci_upper", "n"),
               drop = FALSE]
  class(x_print) <- "data.frame"
  names(x_print)[1] <- mod

  print(x_print, digits = digits, row.names = FALSE)
  invisible(x)
}

#' Compute SE for ATE via influence function
#' @keywords internal
.ate_se <- function(object, cate_values, weights = NULL) {
  n <- length(cate_values)

  if (is.null(weights)) {
    weights <- rep(1, n)
  }
  weights <- weights / sum(weights) * n  # Normalize to sum to n

  ate_est <- .weighted_mean(cate_values, weights)

  # Influence function for ATE: IF_i = tau_i - ATE
  # Under super-population framework
  influence <- (cate_values - ate_est) * weights / n

  # SE = sqrt(sum(IF^2))
  se <- sqrt(sum(influence^2) / n)

  list(
    se = se,
    influence = influence
  )
}
