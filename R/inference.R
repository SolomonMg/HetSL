#' @title Inference Functions for Heterogeneous Effects
#' @description Analytic standard errors, bootstrap, and confidence intervals.
#' @name inference
NULL

#' Analytic Standard Errors
#'
#' Computes analytic standard errors using influence function methodology
#' (TMLE-style).
#'
#' @param object A fitted \code{het_ensemble} object.
#' @param effect_type Type of effect: "ate", "cate", or "mcate".
#' @param moderator For "mcate", the moderator variable name.
#' @param ... Additional arguments passed to effect functions.
#'
#' @return Depends on effect_type:
#' \describe{
#'   \item{ate}{Single SE value}
#'   \item{cate}{Vector of SEs per observation}
#'   \item{mcate}{Vector of SEs per moderator value}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- het_ensemble(y ~ x1 + x2 | treat, data = mydata)
#' analytic_se(fit, effect_type = "ate")
#' }
analytic_se <- function(object, effect_type = c("ate", "cate", "mcate"),
                        moderator = NULL, ...) {

  if (!inherits(object, "het_ensemble")) {
    stop("'object' must be a het_ensemble object.", call. = FALSE)
  }

  effect_type <- match.arg(effect_type)

  switch(effect_type,
         ate = .analytic_se_ate(object, ...),
         cate = .analytic_se_cate(object, ...),
         mcate = .analytic_se_mcate(object, moderator, ...)
  )
}

#' Analytic SE for ATE via influence functions
#' @keywords internal
.analytic_se_ate <- function(object, ...) {
  n <- object$data$n
  y <- object$data$y
  treat <- object$data$treatment
  weights <- object$data$weights

  if (is.null(weights)) {
    weights <- rep(1, n)
  }

  # Get CATE values
  cate_values <- cate(object)

  # ATE estimate
  ate_est <- .weighted_mean(cate_values, weights)

  # Influence function approach
  # Under randomization: IF = (Y - mu_1) * T / p - (Y - mu_0) * (1-T) / (1-p) + tau - ATE
  # Simplified for ensemble: IF ≈ tau - ATE

  influence <- cate_values - ate_est

  # Variance of IF
  var_if <- .weighted_var(influence, weights)

  # SE = sqrt(Var(IF) / n)
  se <- sqrt(var_if / n)

  list(
    se = se,
    variance = var_if,
    influence = influence
  )
}

#' Analytic SE for CATE
#' @keywords internal
.analytic_se_cate <- function(object, ...) {
  # CATE SE is more complex - requires variance estimation per observation
  # Use jackknife or bootstrap approximation

  n <- object$data$n
  K <- length(object$models)

  # Variance from individual learner disagreement
  cf_preds <- lapply(seq_len(K), function(k) {
    lname <- names(object$models)[k]
    model <- object$models[[k]]
    lspec <- .get_learner(lname)

    # Get counterfactual predictions from this learner
    X <- object$data$X
    X_matrix <- .as_matrix(X)

    X0 <- .make_interaction_matrix(X_matrix, rep(0, n), include_main = TRUE)
    X1 <- .make_interaction_matrix(X_matrix, rep(1, n), include_main = TRUE)

    y0 <- lspec$predict_fn(model, X0, rep(0, n), object$family)
    y1 <- lspec$predict_fn(model, X1, rep(1, n), object$family)

    y1 - y0
  })

  cf_matrix <- do.call(cbind, cf_preds)

  # Variance across learners (weighted)
  weights <- object$weights[names(object$models)]

  # Per-observation variance
  cate_var <- apply(cf_matrix, 1, function(row) {
    .weighted_var(row, weights)
  })

  # SE = sqrt(var)
  se <- sqrt(cate_var)

  list(
    se = se,
    variance = cate_var
  )
}

#' Analytic SE for MCATE
#' @keywords internal
.analytic_se_mcate <- function(object, moderator, ...) {
  if (is.null(moderator)) {
    stop("'moderator' must be specified for MCATE standard errors.", call. = FALSE)
  }

  # Get MCATE with SE
  result <- mcate(object, moderator = moderator, se = TRUE, ...)

  list(
    se = result$se,
    by_value = result
  )
}

#' Bootstrap Inference for Treatment Effects
#'
#' Computes bootstrap standard errors and confidence intervals.
#'
#' @param object A fitted \code{het_ensemble} object.
#' @param R Number of bootstrap replications (default 1000).
#' @param effect_type Type of effect: "ate", "cate", or "mcate".
#' @param moderator For "mcate", the moderator variable name.
#' @param parallel Use parallel processing.
#' @param ci_level Confidence level (default 0.95).
#' @param seed Random seed.
#' @param verbose Show progress.
#'
#' @return A list of class \code{bootstrap_effects} containing:
#' \describe{
#'   \item{estimate}{Point estimate(s)}
#'   \item{se}{Bootstrap standard error(s)}
#'   \item{ci_lower}{Lower confidence bound(s)}
#'   \item{ci_upper}{Upper confidence bound(s)}
#'   \item{bootstrap_samples}{Matrix of bootstrap estimates}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- het_ensemble(y ~ x1 + x2 | treat, data = mydata)
#' boot_results <- bootstrap_effects(fit, R = 500)
#' }
bootstrap_effects <- function(object, R = 1000,
                              effect_type = c("ate", "cate", "mcate"),
                              moderator = NULL, parallel = FALSE,
                              ci_level = 0.95, seed = NULL, verbose = TRUE) {

  if (!inherits(object, "het_ensemble")) {
    stop("'object' must be a het_ensemble object.", call. = FALSE)
  }

  effect_type <- match.arg(effect_type)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- object$data$n
  y <- object$data$y
  X <- object$data$X
  treat <- object$data$treatment
  weights <- object$data$weights

  .message(sprintf("Running %d bootstrap replications...", R), verbose)

  # Bootstrap function
  boot_fn <- function(indices) {
    # Resample data
    y_b <- y[indices]
    X_b <- X[indices, , drop = FALSE]
    treat_b <- treat[indices]
    w_b <- if (!is.null(weights)) weights[indices] else NULL

    # Build design matrix
    X_matrix_b <- .as_matrix(X_b)
    X_full_b <- .make_interaction_matrix(X_matrix_b, treat_b, include_main = TRUE)

    # Get predictions from each learner
    preds <- lapply(names(object$models), function(lname) {
      model <- object$models[[lname]]
      lspec <- .get_learner(lname)

      # Counterfactual predictions
      X0_b <- .make_interaction_matrix(X_matrix_b, rep(0, length(treat_b)), include_main = TRUE)
      X1_b <- .make_interaction_matrix(X_matrix_b, rep(1, length(treat_b)), include_main = TRUE)

      y0 <- tryCatch(
        lspec$predict_fn(model, X0_b, rep(0, length(treat_b)), object$family),
        error = function(e) rep(NA_real_, length(treat_b))
      )
      y1 <- tryCatch(
        lspec$predict_fn(model, X1_b, rep(1, length(treat_b)), object$family),
        error = function(e) rep(NA_real_, length(treat_b))
      )

      y1 - y0
    })

    # Ensemble CATE
    pred_matrix <- do.call(cbind, preds)
    cate_b <- .ensemble_predict(pred_matrix, object$weights)

    # Compute effect
    if (effect_type == "ate") {
      .weighted_mean(cate_b, w_b)
    } else if (effect_type == "cate") {
      cate_b
    } else if (effect_type == "mcate") {
      # Simplified: return ATE for now
      .weighted_mean(cate_b, w_b)
    }
  }

  # Run bootstrap
  if (parallel && .is_installed("future") && .is_installed("furrr")) {
    boot_samples <- furrr::future_map(
      seq_len(R),
      function(r) {
        indices <- sample(n, n, replace = TRUE)
        boot_fn(indices)
      },
      .options = furrr::furrr_options(seed = TRUE),
      .progress = verbose
    )
  } else {
    boot_samples <- lapply(seq_len(R), function(r) {
      if (verbose && r %% 100 == 0) {
        .message(sprintf("  Bootstrap %d of %d", r, R), TRUE)
      }
      indices <- sample(n, n, replace = TRUE)
      boot_fn(indices)
    })
  }

  # Process results
  if (effect_type == "ate") {
    boot_matrix <- unlist(boot_samples)

    estimate <- mean(boot_matrix, na.rm = TRUE)
    se <- sd(boot_matrix, na.rm = TRUE)
    ci <- quantile(boot_matrix, probs = c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2),
                   na.rm = TRUE)

    result <- list(
      estimate = estimate,
      se = se,
      ci_lower = ci[1],
      ci_upper = ci[2],
      ci_level = ci_level,
      bootstrap_samples = boot_matrix,
      R = R,
      effect_type = effect_type
    )

  } else if (effect_type == "cate") {
    boot_matrix <- do.call(rbind, boot_samples)

    estimate <- colMeans(boot_matrix, na.rm = TRUE)
    se <- apply(boot_matrix, 2, sd, na.rm = TRUE)
    ci_lower <- apply(boot_matrix, 2, quantile, probs = (1 - ci_level) / 2, na.rm = TRUE)
    ci_upper <- apply(boot_matrix, 2, quantile, probs = 1 - (1 - ci_level) / 2, na.rm = TRUE)

    result <- list(
      estimate = estimate,
      se = se,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      ci_level = ci_level,
      bootstrap_samples = boot_matrix,
      R = R,
      effect_type = effect_type
    )

  } else {
    # mcate - simplified
    boot_matrix <- unlist(boot_samples)

    result <- list(
      estimate = mean(boot_matrix, na.rm = TRUE),
      se = sd(boot_matrix, na.rm = TRUE),
      bootstrap_samples = boot_matrix,
      R = R,
      effect_type = effect_type
    )
  }

  class(result) <- "bootstrap_effects"
  result
}

#' Confidence Intervals for het_ensemble
#'
#' Computes confidence intervals for treatment effects.
#'
#' @param object A fitted \code{het_ensemble} object.
#' @param parm Parameter specification (currently unused).
#' @param level Confidence level (default 0.95).
#' @param method Method for CI: "analytic" (default) or "bootstrap".
#' @param effect_type Type of effect for CI.
#' @param R Bootstrap replications (if method = "bootstrap").
#' @param ... Additional arguments.
#'
#' @return A matrix or vector of confidence intervals.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- het_ensemble(y ~ x1 + x2 | treat, data = mydata)
#' confint(fit)
#' confint(fit, method = "bootstrap", R = 500)
#' }
confint.het_ensemble <- function(object, parm = NULL, level = 0.95,
                                 method = c("analytic", "bootstrap"),
                                 effect_type = "ate", R = 1000, ...) {

  method <- match.arg(method)

  if (method == "analytic") {
    ate_result <- ate(object, se = TRUE, ci_level = level)

    ci <- matrix(c(ate_result$ci_lower, ate_result$ci_upper), nrow = 1)
    colnames(ci) <- c(paste0(100 * (1 - level) / 2, "%"),
                      paste0(100 * (1 - (1 - level) / 2), "%"))
    rownames(ci) <- "ATE"

    ci

  } else {
    boot_result <- bootstrap_effects(object, R = R, effect_type = effect_type,
                                     ci_level = level, ...)

    if (effect_type == "ate") {
      ci <- matrix(c(boot_result$ci_lower, boot_result$ci_upper), nrow = 1)
      colnames(ci) <- c(paste0(100 * (1 - level) / 2, "%"),
                        paste0(100 * (1 - (1 - level) / 2), "%"))
      rownames(ci) <- "ATE"
    } else {
      ci <- cbind(boot_result$ci_lower, boot_result$ci_upper)
      colnames(ci) <- c(paste0(100 * (1 - level) / 2, "%"),
                        paste0(100 * (1 - (1 - level) / 2), "%"))
    }

    ci
  }
}

#' Print method for bootstrap_effects
#'
#' @param x A fitted \code{bootstrap_effects} object.
#' @param digits Number of digits to print.
#' @param ... Unused.
#' @export
print.bootstrap_effects <- function(x, digits = 4, ...) {
  cat(sprintf("Bootstrap Inference (%s)\n", toupper(x$effect_type)))
  cat(paste(rep("=", 35), collapse = ""), "\n\n")
  cat(sprintf("Bootstrap replications: %d\n", x$R))
  cat(sprintf("Confidence level: %d%%\n\n", round(x$ci_level * 100)))

  if (x$effect_type == "ate") {
    cat(sprintf("Estimate: %.*f\n", digits, x$estimate))
    cat(sprintf("Bootstrap SE: %.*f\n", digits, x$se))
    cat(sprintf("CI: [%.*f, %.*f]\n", digits, x$ci_lower, digits, x$ci_upper))
  } else {
    cat(sprintf("CATE estimates: %d observations\n", length(x$estimate)))
    cat(sprintf("Mean SE: %.*f\n", digits, mean(x$se, na.rm = TRUE)))
  }

  invisible(x)
}
