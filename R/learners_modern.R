#' @title Learner Wrappers - Modern HTE Methods
#' @description Fit and predict functions for modern causal ML methods.
#' @name learners_modern
#' @keywords internal
NULL

# =============================================================================
# Causal Forest (grf)
# =============================================================================

#' Fit Causal Forest model
#'
#' Uses the grf package to fit a causal forest.
#' Note: Causal forest directly estimates treatment effects, so we need to
#' convert this to outcome predictions for the ensemble framework.
#'
#' @keywords internal
.fit_causal_forest <- function(X, y, treatment, family = "gaussian",
                               weights = NULL, ...) {
  .require_package("grf", "causal_forest")

  X <- .as_matrix(X)

  # Extract treatment from X if it's included
  # Causal forest needs X without treatment, plus treatment vector

  # Check if 'treat' column exists
  treat_col <- which(colnames(X) == "treat")
  if (length(treat_col) > 0) {
    W <- X[, treat_col]
    # Remove treatment and interaction columns for causal forest
    # Keep only main effect covariates
    interaction_cols <- grep(":treat$", colnames(X))
    remove_cols <- c(treat_col, interaction_cols)
    X_cov <- X[, -remove_cols, drop = FALSE]
  } else {
    W <- treatment
    X_cov <- X
  }

  # Convert treatment to numeric 0/1
  W <- as.numeric(W)

  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }

  # Fit causal forest
  cf <- grf::causal_forest(
    X = X_cov,
    Y = y,
    W = W,
    sample.weights = weights,
    num.trees = 2000,
    honesty = TRUE,
    ...
  )

  # Store additional info
  attr(cf, "family") <- family
  attr(cf, "X_colnames") <- colnames(X_cov)
  attr(cf, "full_X_colnames") <- colnames(X)

  cf
}

#' Predict from Causal Forest model
#'
#' Converts causal forest CATE predictions back to outcome predictions
#' for ensemble compatibility.
#'
#' @keywords internal
.predict_causal_forest <- function(model, newX, treatment = NULL, family = "gaussian") {
  .require_package("grf", "causal_forest")

  newX <- .as_matrix(newX)

  # Extract covariates and treatment from newX
  treat_col <- which(colnames(newX) == "treat")
  if (length(treat_col) > 0) {
    W <- newX[, treat_col]
    interaction_cols <- grep(":treat$", colnames(newX))
    remove_cols <- c(treat_col, interaction_cols)
    X_cov <- newX[, -remove_cols, drop = FALSE]
  } else if (!is.null(treatment)) {
    W <- treatment
    X_cov <- newX
  } else {
    stop("Cannot extract treatment from newX. Provide treatment argument.", call. = FALSE)
  }

  W <- as.numeric(W)

  # Get CATE predictions
  cate_pred <- predict(model, newdata = X_cov)$predictions

  # Get baseline prediction (Y(0) estimate)
  # Use Y.hat from causal forest if available
  if (!is.null(model$Y.hat)) {
    # Y.hat is the estimated E[Y|X] from training
    # For new data, we approximate using the mean
    baseline_pred <- predict(model, newdata = X_cov, estimate.variance = FALSE)

    # Actually, we need to construct outcome predictions
    # Y_pred = E[Y(0)|X] + W * tau(X)
    # We can use: E[Y(0)|X] ≈ E[Y|X] - E[W|X] * tau(X)

    # Simple approximation: use mean outcome adjusted by CATE
    y_mean <- mean(model$Y.orig)
    w_mean <- mean(model$W.orig)

    # E[Y(0)|X] ≈ E[Y|X] - p(W=1) * tau(X)
    # Y_pred(W=w) = E[Y(0)|X] + w * tau(X)
    baseline <- y_mean - w_mean * cate_pred

    y_pred <- baseline + W * cate_pred
  } else {
    # Fallback: simple linear combination
    y_pred <- cate_pred * W
  }

  # For binary outcomes, clip to [0,1]
  if (family == "binomial") {
    y_pred <- .clip(y_pred, 0, 1)
  }

  as.vector(y_pred)
}

# =============================================================================
# XGBoost
# =============================================================================

#' Fit XGBoost model
#' @keywords internal
.fit_xgboost <- function(X, y, treatment, family = "gaussian",
                         weights = NULL, ...) {
  .require_package("xgboost", "xgboost")

  X <- .as_matrix(X)

  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }

  # Determine objective
  if (family == "binomial") {
    objective <- "binary:logistic"
    eval_metric <- "logloss"
  } else {
    objective <- "reg:squarederror"
    eval_metric <- "rmse"
  }

  # Create DMatrix
  dtrain <- xgboost::xgb.DMatrix(
    data = X,
    label = y,
    weight = weights
  )

  # Default parameters
  params <- list(
    objective = objective,
    eval_metric = eval_metric,
    eta = 0.1,
    max_depth = 6,
    subsample = 0.8,
    colsample_bytree = 0.8
  )

  # Merge with user params
  dots <- list(...)
  for (nm in names(dots)) {
    params[[nm]] <- dots[[nm]]
  }

  # Cross-validation to find best nrounds
  cv_result <- xgboost::xgb.cv(
    params = params,
    data = dtrain,
    nrounds = 500,
    nfold = 5,
    early_stopping_rounds = 20,
    verbose = 0
  )

  best_nrounds <- cv_result$best_iteration
  if (is.null(best_nrounds) || best_nrounds < 10) {
    best_nrounds <- 100
  }

  # Train final model
  model <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = best_nrounds,
    verbose = 0
  )

  attr(model, "family") <- family
  attr(model, "colnames") <- colnames(X)

  model
}

#' Predict from XGBoost model
#' @keywords internal
.predict_xgboost <- function(model, newX, treatment = NULL, family = "gaussian") {
  .require_package("xgboost", "xgboost")

  newX <- .as_matrix(newX)

  # Create DMatrix
  dtest <- xgboost::xgb.DMatrix(data = newX)

  # Predict
  preds <- predict(model, dtest)

  # For binary:logistic, predictions are already probabilities
  as.vector(preds)
}

# =============================================================================
# Helper for causal forest integration
# =============================================================================

#' Get CATE directly from causal forest
#'
#' This is useful for treatment effect computation, bypassing the
#' ensemble prediction conversion.
#'
#' @param model Causal forest model
#' @param newX New covariate matrix
#' @return Vector of CATE estimates
#' @keywords internal
.get_causal_forest_cate <- function(model, newX) {
  .require_package("grf", "causal_forest")

  newX <- .as_matrix(newX)

  # Extract just the covariates
  treat_col <- which(colnames(newX) == "treat")
  if (length(treat_col) > 0) {
    interaction_cols <- grep(":treat$", colnames(newX))
    remove_cols <- c(treat_col, interaction_cols)
    X_cov <- newX[, -remove_cols, drop = FALSE]
  } else {
    X_cov <- newX
  }

  pred <- predict(model, newdata = X_cov, estimate.variance = TRUE)

  list(
    cate = pred$predictions,
    variance = pred$variance.estimates
  )
}

#' Get CATE variance from causal forest
#' @keywords internal
.get_causal_forest_variance <- function(model, newX) {
  result <- .get_causal_forest_cate(model, newX)
  result$variance
}
