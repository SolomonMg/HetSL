#' @title Prediction Methods for het_ensemble
#' @description Predict outcomes from fitted ensemble models.
#' @name predict.het_ensemble
NULL

#' Predict from Heterogeneous Effects Ensemble
#'
#' Generate predictions from a fitted \code{het_ensemble} model.
#'
#' @param object A fitted \code{het_ensemble} object.
#' @param newdata Optional data frame of new observations. If NULL, returns
#'   predictions for training data.
#' @param type Type of prediction:
#'   \describe{
#'     \item{"ensemble"}{Weighted ensemble prediction (default)}
#'     \item{"individual"}{Predictions from each learner}
#'     \item{"cate"}{Conditional average treatment effect}
#'     \item{"counterfactual"}{Predictions under both treatment conditions}
#'   }
#' @param treatment For counterfactual predictions, the treatment value(s) to use.
#'   If NULL, uses actual treatment values in data.
#' @param ... Additional arguments (unused).
#'
#' @return Depends on \code{type}:
#'   \describe{
#'     \item{"ensemble"}{Numeric vector of predictions}
#'     \item{"individual"}{Matrix with columns for each learner}
#'     \item{"cate"}{Numeric vector of treatment effect estimates}
#'     \item{"counterfactual"}{Data frame with columns: y0, y1, cate}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Ensemble predictions
#' pred <- predict(fit)
#'
#' # Individual learner predictions
#' pred_all <- predict(fit, type = "individual")
#'
#' # CATE estimates
#' cate <- predict(fit, type = "cate")
#'
#' # Counterfactual predictions
#' cf <- predict(fit, type = "counterfactual")
#' }
predict.het_ensemble <- function(object, newdata = NULL,
                                 type = c("ensemble", "individual", "cate", "counterfactual"),
                                 treatment = NULL, ...) {

  type <- match.arg(type)

  # Use training data if newdata not provided
  if (is.null(newdata)) {
    X <- object$data$X
    treat <- object$data$treatment
    X_full <- object$data$X_full
    use_training <- TRUE
  } else {
    # Build design matrix for new data
    X <- newdata[, object$covariate_vars, drop = FALSE]

    # Get treatment values
    if (!is.null(treatment)) {
      if (is.matrix(treatment) || is.data.frame(treatment)) {
        treat <- treatment
      } else {
        treat <- rep(treatment, nrow(X))[seq_len(nrow(X))]
      }
    } else if (all(object$treatment_var %in% names(newdata))) {
      treat <- newdata[, object$treatment_var, drop = FALSE]
      treat <- .recode_treatment(treat)
    } else {
      # No treatment - will need to compute counterfactuals
      treat <- NULL
    }

    use_training <- FALSE
    X_full <- NULL
  }

  # Route to appropriate prediction function
  switch(type,
         ensemble = .predict_ensemble(object, X, treat, X_full),
         individual = .predict_individual(object, X, treat, X_full),
         cate = .predict_cate(object, X, X_full),
         counterfactual = .predict_counterfactual(object, X, X_full, treatment)
  )
}

#' Ensemble prediction
#' @keywords internal
.predict_ensemble <- function(object, X, treatment, X_full = NULL) {
  # Get individual predictions
  pred_matrix <- .predict_individual(object, X, treatment, X_full)

  # Compute weighted average
  .ensemble_predict(pred_matrix, object$weights)
}

#' Individual learner predictions
#' @keywords internal
.predict_individual <- function(object, X, treatment, X_full = NULL) {
  # Build design matrix if needed
  if (is.null(X_full)) {
    if (is.null(treatment)) {
      stop("Treatment must be provided for prediction.", call. = FALSE)
    }
    X_matrix <- .as_matrix(X)
    X_full <- .make_interaction_matrix(X_matrix, treatment, include_main = TRUE)
  }

  n <- nrow(X_full)
  n_learners <- length(object$models)

  # Get predictions from each learner
  pred_matrix <- matrix(NA_real_, nrow = n, ncol = n_learners)
  colnames(pred_matrix) <- names(object$models)

  for (i in seq_along(object$models)) {
    lname <- names(object$models)[i]
    model <- object$models[[i]]
    lspec <- .get_learner(lname)

    preds <- tryCatch(
      lspec$predict_fn(model, X_full, treatment, object$family),
      error = function(e) {
        warning(sprintf("Prediction failed for %s: %s", lname, e$message),
                call. = FALSE)
        rep(NA_real_, n)
      }
    )

    pred_matrix[, i] <- preds
  }

  pred_matrix
}

#' CATE prediction (treatment effect for each observation)
#' @keywords internal
.predict_cate <- function(object, X, X_full = NULL) {
  if (is.matrix(object$data$treatment)) {
    stop("CATE is only supported for binary treatments.", call. = FALSE)
  }
  # Compute counterfactual predictions
  cf <- .predict_counterfactual(object, X, X_full)
  cf$cate
}

#' Counterfactual predictions
#' @keywords internal
.predict_counterfactual <- function(object, X, X_full = NULL, treatment = NULL) {
  if (is.matrix(object$data$treatment)) {
    if (is.list(treatment) &&
        !is.null(treatment$treated) &&
        !is.null(treatment$control)) {
      X_matrix <- .as_matrix(X)
      pred_control <- .predict_individual(object, X, treatment$control, X_full)
      y0 <- .ensemble_predict(pred_control, object$weights)

      pred_treated <- .predict_individual(object, X, treatment$treated, X_full)
      y1 <- .ensemble_predict(pred_treated, object$weights)

      return(data.frame(y0 = y0, y1 = y1, cate = y1 - y0))
    }

    stop("Counterfactual predictions require treatment=list(control=..., treated=...) for multi-treatment models.",
         call. = FALSE)
  }

  X_matrix <- .as_matrix(X)
  n <- nrow(X_matrix)

  # Prediction under control (treatment = 0)
  X_control <- .make_interaction_matrix(X_matrix, rep(0, n), include_main = TRUE)
  pred_control <- .predict_individual(object, X, rep(0, n), X_control)
  y0 <- .ensemble_predict(pred_control, object$weights)

  # Prediction under treatment (treatment = 1)
  X_treated <- .make_interaction_matrix(X_matrix, rep(1, n), include_main = TRUE)
  pred_treated <- .predict_individual(object, X, rep(1, n), X_treated)
  y1 <- .ensemble_predict(pred_treated, object$weights)

  # CATE = E[Y(1) - Y(0) | X]
  cate <- y1 - y0

  data.frame(
    y0 = y0,
    y1 = y1,
    cate = cate
  )
}

#' Get fitted values
#'
#' @param object het_ensemble object
#' @param type Type of fitted values
#' @return Vector of fitted values
#' @keywords internal
fitted.het_ensemble <- function(object, type = c("ensemble", "individual"), ...) {
  type <- match.arg(type)

  if (type == "ensemble") {
    .ensemble_predict(object$cv_predictions, object$weights)
  } else {
    object$cv_predictions
  }
}

#' Get residuals
#'
#' @param object het_ensemble object
#' @param type Type of residuals
#' @return Vector of residuals
#' @keywords internal
residuals.het_ensemble <- function(object, type = c("response", "deviance"), ...) {
  type <- match.arg(type)

  y <- object$data$y
  fitted <- fitted.het_ensemble(object, type = "ensemble")

  if (type == "response") {
    y - fitted
  } else {
    # Deviance residuals
    if (object$family == "binomial") {
      fitted <- .clip(fitted, 1e-10, 1 - 1e-10)
      sign(y - fitted) * sqrt(-2 * (y * log(fitted) + (1 - y) * log(1 - fitted)))
    } else {
      y - fitted
    }
  }
}

#' Batch prediction for large datasets
#'
#' Predicts in batches to manage memory for large datasets.
#'
#' @param object het_ensemble object
#' @param newdata Data frame
#' @param batch_size Number of observations per batch
#' @param type Prediction type
#' @param verbose Show progress
#' @return Predictions
#' @keywords internal
.predict_batch <- function(object, newdata, batch_size = 10000,
                           type = "ensemble", verbose = TRUE) {
  n <- nrow(newdata)

  if (n <= batch_size) {
    return(predict(object, newdata = newdata, type = type))
  }

  # Split into batches
  n_batches <- ceiling(n / batch_size)
  batch_starts <- seq(1, n, by = batch_size)
  batch_ends <- pmin(batch_starts + batch_size - 1, n)

  results <- vector("list", n_batches)

  for (i in seq_len(n_batches)) {
    if (verbose) {
      .message(sprintf("Batch %d of %d...", i, n_batches), TRUE)
    }

    batch_data <- newdata[batch_starts[i]:batch_ends[i], , drop = FALSE]
    results[[i]] <- predict(object, newdata = batch_data, type = type)
  }

  # Combine results
  if (type %in% c("ensemble", "cate")) {
    unlist(results)
  } else if (type == "individual") {
    do.call(rbind, results)
  } else if (type == "counterfactual") {
    do.call(rbind, results)
  } else {
    results
  }
}
