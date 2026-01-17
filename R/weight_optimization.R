#' @title Weight Optimization for Ensemble
#' @description Quadratic programming for optimal ensemble weights.
#' @name weight_optimization
#' @keywords internal
NULL

#' Optimize ensemble weights via quadratic programming
#'
#' Finds optimal convex combination weights that minimize CV loss.
#' Weights are constrained to be non-negative and sum to 1.
#'
#' @param cv_predictions Matrix of CV predictions (n x K)
#' @param y True outcomes
#' @param family "binomial" or "gaussian"
#' @param weights Observation weights
#' @param lower_bound Lower bound for weights (default 0)
#' @return Named vector of optimal weights
#' @keywords internal
#'
#' @details
#' Solves the quadratic program:
#' \deqn{\min_\alpha \frac{1}{2} \alpha' D \alpha - d' \alpha}
#' subject to:
#' \deqn{\sum_k \alpha_k = 1}
#' \deqn{\alpha_k \geq 0}
#'
#' Uses the \code{quadprog} package for the optimization.
.optimize_weights <- function(cv_predictions, y, family = "gaussian",
                              weights = NULL, lower_bound = 0) {

  K <- ncol(cv_predictions)
  learner_names <- colnames(cv_predictions)

  if (K == 1) {
    # Single learner - weight is 1
    result <- setNames(1.0, learner_names)
    return(result)
  }

  # Get QP components
  qp_components <- .cv_risk_matrix(cv_predictions, y, family, weights)
  Dmat <- qp_components$Dmat
  dvec <- qp_components$dvec

  # Constraints:
  # 1. Sum to 1: sum(alpha) = 1
  # 2. Non-negativity: alpha_k >= lower_bound

  # Amat is (K x m) where m is number of constraints
  # First column: equality constraint (sum = 1)
  # Remaining columns: inequality constraints (alpha_k >= lower_bound)

  # quadprog solves: min 0.5 * x' D x - d' x
  # subject to: A' x >= b
  # First meq constraints are equalities

  # Equality: sum(alpha) = 1 -> [1,1,...,1]' alpha = 1
  # Inequality: alpha_k >= 0 -> I alpha >= 0

  Amat <- cbind(
    rep(1, K),           # Sum to 1 constraint
    diag(1, K)           # Non-negativity constraints
  )

  bvec <- c(
    1,                   # Sum equals 1
    rep(lower_bound, K)  # Lower bounds
  )

  # Solve QP
  qp_result <- tryCatch(
    quadprog::solve.QP(
      Dmat = Dmat,
      dvec = dvec,
      Amat = Amat,
      bvec = bvec,
      meq = 1  # First constraint is equality
    ),
    error = function(e) {
      # Try with slightly more regularization
      Dmat_reg <- Dmat + diag(1e-6, K)
      quadprog::solve.QP(
        Dmat = Dmat_reg,
        dvec = dvec,
        Amat = Amat,
        bvec = bvec,
        meq = 1
      )
    }
  )

  # Extract weights
  alpha <- qp_result$solution

  # Clean up numerical issues
  alpha[alpha < 1e-10] <- 0
  alpha <- alpha / sum(alpha)  # Ensure exactly sums to 1

  names(alpha) <- learner_names
  alpha
}

#' Alternative weight optimization via NNLS
#'
#' Non-negative least squares approach when QP fails.
#'
#' @param cv_predictions Matrix of CV predictions
#' @param y True outcomes
#' @param family "binomial" or "gaussian"
#' @param weights Observation weights
#' @return Named vector of weights
#' @keywords internal
.optimize_weights_nnls <- function(cv_predictions, y, family = "gaussian",
                                   weights = NULL) {
  K <- ncol(cv_predictions)
  learner_names <- colnames(cv_predictions)

  if (is.null(weights)) {
    weights <- rep(1, nrow(cv_predictions))
  }

  W <- sqrt(weights)
  X <- sweep(cv_predictions, 1, W, "*")
  y_w <- y * W

  # Solve using coordinate descent
  alpha <- rep(1/K, K)
  max_iter <- 1000
  tol <- 1e-8

  for (iter in seq_len(max_iter)) {
    alpha_old <- alpha

    for (k in seq_len(K)) {
      # Residual without contribution from k
      resid <- y_w - X[, -k, drop = FALSE] %*% alpha[-k]

      # Optimal alpha_k (unconstrained)
      alpha_k <- sum(X[, k] * resid) / sum(X[, k]^2)

      # Project to non-negative
      alpha[k] <- max(0, alpha_k)
    }

    # Check convergence
    if (max(abs(alpha - alpha_old)) < tol) {
      break
    }
  }

  # Normalize
  if (sum(alpha) > 0) {
    alpha <- alpha / sum(alpha)
  } else {
    alpha <- rep(1/K, K)
  }

  names(alpha) <- learner_names
  alpha
}

#' Compute ensemble predictions given weights
#'
#' @param predictions Matrix of predictions (n x K) or list of prediction vectors
#' @param weights Named vector of weights
#' @return Vector of ensemble predictions
#' @keywords internal
.ensemble_predict <- function(predictions, weights) {
  if (is.list(predictions)) {
    predictions <- do.call(cbind, predictions)
  }

  # Align weights with prediction columns
  if (!is.null(names(weights)) && !is.null(colnames(predictions))) {
    common <- intersect(names(weights), colnames(predictions))
    weights <- weights[common]
    predictions <- predictions[, common, drop = FALSE]
  }

  # Normalize weights
  weights <- weights / sum(weights)

  # Compute weighted average
  as.vector(predictions %*% weights)
}

#' Get CV-based performance metrics for learners
#'
#' @param cv_predictions Matrix of CV predictions
#' @param y True outcomes
#' @param family "binomial" or "gaussian"
#' @param weights Observation weights
#' @return Data frame with learner performance
#' @keywords internal
.learner_performance <- function(cv_predictions, y, family = "gaussian",
                                 weights = NULL) {
  K <- ncol(cv_predictions)
  learner_names <- colnames(cv_predictions)
  n <- length(y)

  if (is.null(weights)) {
    weights <- rep(1, n)
  }

  metrics <- data.frame(
    learner = learner_names,
    stringsAsFactors = FALSE
  )

  if (family == "binomial") {
    # Log loss
    metrics$cv_logloss <- apply(cv_predictions, 2, function(preds) {
      preds <- .clip(preds, 1e-10, 1 - 1e-10)
      -sum(weights * (y * log(preds) + (1 - y) * log(1 - preds))) / sum(weights)
    })

    # Accuracy
    metrics$cv_accuracy <- apply(cv_predictions, 2, function(preds) {
      pred_class <- as.numeric(preds > 0.5)
      sum(weights * (pred_class == y)) / sum(weights)
    })

    # AUC (simple version)
    metrics$cv_auc <- apply(cv_predictions, 2, function(preds) {
      .simple_auc(y, preds, weights)
    })

  } else {
    # MSE
    metrics$cv_mse <- apply(cv_predictions, 2, function(preds) {
      sum(weights * (y - preds)^2) / sum(weights)
    })

    # RMSE
    metrics$cv_rmse <- sqrt(metrics$cv_mse)

    # MAE
    metrics$cv_mae <- apply(cv_predictions, 2, function(preds) {
      sum(weights * abs(y - preds)) / sum(weights)
    })

    # R-squared
    ss_tot <- sum(weights * (y - .weighted_mean(y, weights))^2)
    metrics$cv_r2 <- 1 - metrics$cv_mse * sum(weights) / ss_tot
  }

  metrics
}

#' Simple AUC calculation
#' @keywords internal
.simple_auc <- function(y, pred, weights = NULL) {
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }

  # Sort by predictions
  ord <- order(pred, decreasing = TRUE)
  y <- y[ord]
  weights <- weights[ord]

  # Weighted AUC
  n1 <- sum(weights[y == 1])
  n0 <- sum(weights[y == 0])

  if (n1 == 0 || n0 == 0) return(0.5)

  # Cumulative sum of weights for class 0
  cum_w0 <- cumsum(weights * (1 - y))

  # AUC = sum of weights of 1s * cum_w0 at each 1
  auc <- sum(weights[y == 1] * cum_w0[y == 1]) / (n1 * n0)

  # This gives ROC-AUC
  1 - auc
}
