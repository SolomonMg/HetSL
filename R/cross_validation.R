#' @title Cross-Validation for Ensemble Learning
#' @description D-fold cross-validation infrastructure for learner evaluation.
#' @name cross_validation
#' @keywords internal
NULL

#' Create stratified cross-validation folds
#'
#' Creates fold assignments stratified by treatment to ensure balance.
#'
#' @param n Number of observations
#' @param nfolds Number of folds (default 10)
#' @param treatment Treatment vector for stratification (optional)
#' @param seed Random seed
#' @return Integer vector of fold assignments
#' @keywords internal
.create_folds <- function(n, nfolds = 10L, treatment = NULL, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.matrix(treatment) || is.data.frame(treatment)) {
    treatment <- NULL
  }

  if (is.null(treatment)) {
    # Simple random assignment
    fold_ids <- sample(rep(seq_len(nfolds), length.out = n))
  } else {
    # Stratified by treatment
    fold_ids <- integer(n)
    treatment <- as.factor(treatment)

    for (level in levels(treatment)) {
      idx <- which(treatment == level)
      n_level <- length(idx)
      fold_ids[idx] <- sample(rep(seq_len(nfolds), length.out = n_level))
    }
  }

  fold_ids
}

#' Run cross-validation for a single learner
#'
#' @param X Covariate matrix (with interactions)
#' @param y Outcome vector
#' @param treatment Treatment vector
#' @param folds Fold assignments
#' @param learner_name Name of learner
#' @param fit_fn Fitting function
#' @param predict_fn Prediction function
#' @param family "binomial" or "gaussian"
#' @param weights Observation weights
#' @param on_error Error handling mode
#' @param verbose Print progress
#' @param ... Additional arguments to fit function
#' @return List with cv_predictions and learner_info
#' @keywords internal
.cv_single_learner <- function(X, y, treatment, folds, learner_name,
                               fit_fn, predict_fn, family = "gaussian",
                               weights = NULL, on_error = "warn",
                               verbose = TRUE, ...) {

  n <- length(y)
  nfolds <- max(folds)
  cv_predictions <- rep(NA_real_, n)
  fold_models <- vector("list", nfolds)
  success <- TRUE
  error_msg <- NULL

  .message(sprintf("  Fitting %s...", learner_name), verbose)

  for (fold in seq_len(nfolds)) {
    train_idx <- which(folds != fold)
    test_idx <- which(folds == fold)

    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    treat_train <- treatment[train_idx]
    X_test <- X[test_idx, , drop = FALSE]
    treat_test <- treatment[test_idx]

    w_train <- if (!is.null(weights)) weights[train_idx] else NULL

    # Fit model
    model <- .safe_try(
      fit_fn(X_train, y_train, treat_train, family = family,
             weights = w_train, ...),
      on_error = on_error,
      learner_name = sprintf("%s (fold %d)", learner_name, fold)
    )

    if (is.null(model)) {
      success <- FALSE
      error_msg <- sprintf("Failed to fit %s in fold %d", learner_name, fold)
      break
    }

    # Predict on held-out fold
    preds <- .safe_try(
      predict_fn(model, X_test, treat_test, family = family),
      on_error = on_error,
      learner_name = sprintf("%s predict (fold %d)", learner_name, fold)
    )

    if (is.null(preds)) {
      success <- FALSE
      error_msg <- sprintf("Failed to predict with %s in fold %d", learner_name, fold)
      break
    }

    cv_predictions[test_idx] <- preds
    fold_models[[fold]] <- model
  }

  if (!success) {
    return(list(
      cv_predictions = NULL,
      success = FALSE,
      error = error_msg
    ))
  }

  list(
    cv_predictions = cv_predictions,
    fold_models = fold_models,
    learner_name = learner_name,
    success = TRUE,
    error = NULL
  )
}

#' Run cross-validation for all learners
#'
#' @param X Covariate matrix
#' @param y Outcome vector
#' @param treatment Treatment vector
#' @param learners Named list of learner specs
#' @param nfolds Number of folds
#' @param family "binomial" or "gaussian"
#' @param weights Observation weights
#' @param on_error Error handling
#' @param parallel Use parallel processing
#' @param seed Random seed
#' @param verbose Print progress
#' @return List with cv_predictions matrix and learner results
#' @keywords internal
.run_cv_all_learners <- function(X, y, treatment, learners, nfolds = 10L,
                                 family = "gaussian", weights = NULL,
                                 on_error = "warn", parallel = FALSE,
                                 seed = NULL, verbose = TRUE) {

  n <- length(y)
  n_learners <- length(learners)

  if (nfolds <= 1) {
    .message("Running single-fold (no CV) fit for all learners...", verbose)

    fit_one <- function(lname) {
      lspec <- learners[[lname]]
      fit_time <- system.time({
        model <- .safe_try(
          lspec$fit_fn(X, y, treatment, family = family, weights = weights),
          on_error = on_error,
          learner_name = lname
        )
      })

      if (is.null(model)) {
        return(list(success = FALSE, cv_predictions = NULL, fit_time = fit_time))
      }

      pred_time <- system.time({
        preds <- .safe_try(
          lspec$predict_fn(model, X, treatment, family = family),
          on_error = on_error,
          learner_name = sprintf("%s predict", lname)
        )
      })

      if (is.null(preds)) {
        return(list(success = FALSE, cv_predictions = NULL, fit_time = fit_time))
      }

      list(
        success = TRUE,
        cv_predictions = preds,
        fit_time = fit_time,
        pred_time = pred_time
      )
    }

    cv_results <- lapply(names(learners), fit_one)
    names(cv_results) <- names(learners)

    successful <- sapply(cv_results, function(x) x$success)
    if (sum(successful) == 0) {
      stop("All learners failed during single-fold fit.", call. = FALSE)
    }

    successful_names <- names(cv_results)[successful]
    cv_pred_matrix <- do.call(cbind, lapply(successful_names, function(lname) {
      cv_results[[lname]]$cv_predictions
    }))
    colnames(cv_pred_matrix) <- successful_names

    folds <- rep(1L, n)
    return(list(
      cv_predictions = cv_pred_matrix,
      cv_results = cv_results[successful],
      folds = folds,
      successful_learners = successful_names,
      failed_learners = names(cv_results)[!successful]
    ))
  }

  # Create folds
  folds <- .create_folds(n, nfolds, treatment, seed)

  .message(sprintf("Running %d-fold CV with %d learners...", nfolds, n_learners),
           verbose)

  # Run CV for each learner
  if (parallel && .is_installed("future") && .is_installed("furrr")) {
    # Parallel execution
    cv_results <- furrr::future_map(
      names(learners),
      function(lname) {
        lspec <- learners[[lname]]
        .cv_single_learner(
          X = X, y = y, treatment = treatment, folds = folds,
          learner_name = lname,
          fit_fn = lspec$fit_fn,
          predict_fn = lspec$predict_fn,
          family = family,
          weights = weights,
          on_error = on_error,
          verbose = FALSE
        )
      },
      .options = furrr::furrr_options(seed = TRUE)
    )
    names(cv_results) <- names(learners)

  } else {
    # Sequential execution
    cv_results <- lapply(names(learners), function(lname) {
      lspec <- learners[[lname]]
      .cv_single_learner(
        X = X, y = y, treatment = treatment, folds = folds,
        learner_name = lname,
        fit_fn = lspec$fit_fn,
        predict_fn = lspec$predict_fn,
        family = family,
        weights = weights,
        on_error = on_error,
        verbose = verbose
      )
    })
    names(cv_results) <- names(learners)
  }

  # Filter successful learners
  successful <- sapply(cv_results, function(x) x$success)

  if (sum(successful) == 0) {
    stop("All learners failed during cross-validation.", call. = FALSE)
  }

  if (sum(!successful) > 0 && on_error != "ignore") {
    failed_names <- names(cv_results)[!successful]
    warning(sprintf("Learners failed and were excluded: %s",
                    paste(failed_names, collapse = ", ")), call. = FALSE)
  }

  # Build prediction matrix from successful learners
  successful_names <- names(cv_results)[successful]
  cv_pred_matrix <- do.call(cbind, lapply(successful_names, function(lname) {
    cv_results[[lname]]$cv_predictions
  }))
  colnames(cv_pred_matrix) <- successful_names

  list(
    cv_predictions = cv_pred_matrix,
    cv_results = cv_results[successful],
    folds = folds,
    successful_learners = successful_names,
    failed_learners = names(cv_results)[!successful]
  )
}

#' Calculate CV loss for a prediction matrix
#'
#' @param cv_predictions Matrix of CV predictions (n x n_learners)
#' @param y True outcomes
#' @param family "binomial" or "gaussian"
#' @param weights Observation weights
#' @return Vector of CV losses by learner
#' @keywords internal
.cv_loss <- function(cv_predictions, y, family = "gaussian", weights = NULL) {
  n <- length(y)
  if (is.null(weights)) {
    weights <- rep(1, n)
  }

  if (family == "binomial") {
    # Log loss (cross-entropy)
    apply(cv_predictions, 2, function(preds) {
      preds <- .clip(preds, 1e-10, 1 - 1e-10)
      -sum(weights * (y * log(preds) + (1 - y) * log(1 - preds))) / sum(weights)
    })
  } else {
    # MSE
    apply(cv_predictions, 2, function(preds) {
      sum(weights * (y - preds)^2) / sum(weights)
    })
  }
}

#' Calculate CV risk matrix for ensemble weights
#'
#' Computes the quadratic programming components for weight optimization.
#'
#' @param cv_predictions Matrix of CV predictions
#' @param y True outcomes
#' @param family "binomial" or "gaussian"
#' @param weights Observation weights
#' @return List with Dmat and dvec for QP
#' @keywords internal
.cv_risk_matrix <- function(cv_predictions, y, family = "gaussian",
                            weights = NULL) {
  n <- nrow(cv_predictions)
  K <- ncol(cv_predictions)

  if (is.null(weights)) {
    weights <- rep(1, n)
  }

  if (family == "binomial") {
    # Transform to working residuals for log-loss minimization
    # Use pseudo-residuals approach
    cv_predictions <- .clip(cv_predictions, 1e-10, 1 - 1e-10)

    # Gradient of log-loss
    residuals <- y - cv_predictions

    # Dmat: cross-product of weighted residuals
    # This approximates the Hessian
    W <- sqrt(weights)
    weighted_resid <- sweep(residuals, 1, W, "*")
    Dmat <- crossprod(weighted_resid)

    # dvec: weighted sum of residuals * predictions
    dvec <- colSums(sweep(residuals * cv_predictions, 1, weights, "*"))

  } else {
    # Squared error loss
    # Objective: min_w sum_i w_i (y_i - sum_k alpha_k pred_ik)^2
    # = alpha' Dmat alpha - 2 dvec' alpha + const

    W <- sqrt(weights)
    weighted_preds <- sweep(cv_predictions, 1, W, "*")
    Dmat <- crossprod(weighted_preds)

    weighted_y <- y * W
    dvec <- crossprod(weighted_preds, weighted_y)
  }

  # Add small ridge for numerical stability
  Dmat <- Dmat + diag(1e-8, K)

  list(
    Dmat = Dmat,
    dvec = as.vector(dvec)
  )
}
