#' @title Learner Wrappers - Original Paper Methods
#' @description Fit and predict functions for learners from Grimmer, Messing & Westwood (2017).
#' @name learners
#' @keywords internal
NULL

# =============================================================================
# LASSO (glmnet)
# =============================================================================

#' Fit LASSO model
#' @keywords internal
.fit_lasso <- function(X, y, treatment, family = "gaussian", weights = NULL, ...) {
  X <- .as_matrix(X)

  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }

  # Use cv.glmnet to select lambda
  fam <- if (family == "binomial") "binomial" else "gaussian"

  model <- glmnet::cv.glmnet(
    x = X,
    y = y,
    family = fam,
    alpha = 1,  # LASSO
    weights = weights,
    nfolds = 5,
    ...
  )

  model
}

#' Predict from LASSO model
#' @keywords internal
.predict_lasso <- function(model, newX, treatment = NULL, family = "gaussian") {
  newX <- .as_matrix(newX)

  type <- if (family == "binomial") "response" else "response"

  as.vector(predict(model, newx = newX, s = "lambda.min", type = type))
}

# =============================================================================
# Elastic Net (glmnet)
# =============================================================================

#' Fit Elastic Net model
#' @keywords internal
.fit_elastic_net <- function(X, y, treatment, family = "gaussian",
                             weights = NULL, alpha = 0.5, ...) {
  X <- .as_matrix(X)

  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }

  fam <- if (family == "binomial") "binomial" else "gaussian"

  model <- glmnet::cv.glmnet(
    x = X,
    y = y,
    family = fam,
    alpha = alpha,
    weights = weights,
    nfolds = 5,
    ...
  )

  attr(model, "alpha") <- alpha
  model
}

#' Predict from Elastic Net model
#' @keywords internal
.predict_elastic_net <- function(model, newX, treatment = NULL, family = "gaussian") {
  newX <- .as_matrix(newX)

  type <- if (family == "binomial") "response" else "response"

  as.vector(predict(model, newx = newX, s = "lambda.min", type = type))
}

# =============================================================================
# BayesGLM (arm)
# =============================================================================

#' Fit Bayesian GLM
#' @keywords internal
.fit_bayesglm <- function(X, y, treatment, family = "gaussian", weights = NULL, ...) {
  .require_package("arm", "bayesglm")

  X <- .as_matrix(X)

  # Create data frame for arm::bayesglm
  df <- as.data.frame(X)
  df$.y <- y

  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }

  # Formula
  covar_names <- colnames(X)
  # Sanitize names for formula
  safe_names <- make.names(covar_names, unique = TRUE)
  colnames(df)[seq_along(covar_names)] <- safe_names

  fml <- as.formula(paste(".y ~", paste(safe_names, collapse = " + ")))

  fam <- if (family == "binomial") binomial() else gaussian()

  model <- arm::bayesglm(
    formula = fml,
    family = fam,
    data = df,
    weights = weights,
    ...
  )

  # Store name mapping
  attr(model, "covar_mapping") <- setNames(safe_names, covar_names)

  model
}

#' Predict from Bayesian GLM
#' @keywords internal
.predict_bayesglm <- function(model, newX, treatment = NULL, family = "gaussian") {
  .require_package("arm", "bayesglm")

  newX <- .as_matrix(newX)

  # Map column names
  mapping <- attr(model, "covar_mapping")
  if (!is.null(mapping)) {
    orig_names <- names(mapping)
    safe_names <- unname(mapping)
    # Try to match by position if names don't match
    if (!is.null(colnames(newX)) && all(colnames(newX) %in% orig_names)) {
      colnames(newX) <- mapping[colnames(newX)]
    } else if (ncol(newX) == length(safe_names)) {
      colnames(newX) <- safe_names
    }
  }

  newdf <- as.data.frame(newX)

  type <- if (family == "binomial") "response" else "response"

  as.vector(predict(model, newdata = newdf, type = type))
}

# =============================================================================
# FindIt
# =============================================================================

#' Fit FindIt model
#' @keywords internal
.fit_findit <- function(X, y, treatment, family = "binomial", weights = NULL, ...) {
  .require_package("FindIt", "FindIt")

  if (family != "binomial") {
    stop("FindIt is only supported for binomial outcomes.", call. = FALSE)
  }

  X <- .as_matrix(X)

  # Extract covariates and treatment components from interaction matrix.
  treat_col <- which(colnames(X) == "treat")
  interaction_cols <- grep(":treat$", colnames(X))

  if (!is.null(treatment)) {
    X_cov <- X[, setdiff(seq_len(ncol(X)), c(treat_col, interaction_cols)), drop = FALSE]
    treat_mat <- as.matrix(treatment)
  } else if (length(treat_col) > 0) {
    X_cov <- X[, setdiff(seq_len(ncol(X)), c(treat_col, interaction_cols)), drop = FALSE]
    treat_mat <- as.matrix(X[, treat_col, drop = FALSE])
  } else {
    stop("FindIt requires treatment information in X or via 'treatment'.", call. = FALSE)
  }

  std <- .standardize_matrix(X_cov, center = TRUE, scale = TRUE)
  Xstd <- std$X

  # FindIt expects Y in {-1, 1}.
  FIY <- y
  FIY[FIY == 0] <- -1

  if (is.null(colnames(Xstd))) {
    colnames(Xstd) <- paste0("Cov", seq_len(ncol(Xstd)))
  }
  if (is.null(colnames(treat_mat))) {
    colnames(treat_mat) <- paste0("treat", seq_len(ncol(treat_mat)))
  }

  treat_mat <- as.matrix(treat_mat)

  start <- model.matrix(~Xstd * treat_mat)
  treat_start <- grep("^treat", colnames(start))[1]
  if (is.na(treat_start)) {
    stop("FindIt could not identify treatment columns.", call. = FALSE)
  }
  treat2 <- start[, treat_start:ncol(start), drop = FALSE]

  call_findit <- function(use_xc = TRUE) {
    x_arg <- if (use_xc) list(X.c = Xstd) else list(X = Xstd)
    if (ncol(treat_mat) > 1) {
      do.call(
        FindIt::FindIt,
        c(
          list(FIY),
          x_arg,
          list(
            treat2,
            type = "multiple",
            scale.c = std$scale,
            search.lambdas = TRUE,
            fit.glmnet = TRUE,
            wts = if (is.null(weights)) 1 else weights
          ),
          list(...)
        )
      )
    } else {
      do.call(
        FindIt::FindIt,
        c(
          list(FIY),
          x_arg,
          list(
            treat_mat,
            type = "single",
            scale.c = std$scale,
            search.lambdas = TRUE,
            fit.glmnet = TRUE,
            wts = if (is.null(weights)) 1 else weights
          ),
          list(...)
        )
      )
    }
  }

  fit <- tryCatch(
    call_findit(TRUE),
    error = function(e) call_findit(FALSE)
  )

  attr(fit, "findit_center") <- std$center
  attr(fit, "findit_scale") <- std$scale
  attr(fit, "findit_treat_start") <- treat_start
  attr(fit, "findit_treat_cols") <- colnames(treat2)
  attr(fit, "family") <- family
  fit
}

#' Predict from FindIt model
#' @keywords internal
.predict_findit <- function(model, newX, treatment = NULL, family = "binomial") {
  .require_package("FindIt", "FindIt")

  newX <- .as_matrix(newX)

  treat_col <- which(colnames(newX) == "treat")
  interaction_cols <- grep(":treat$", colnames(newX))

  if (!is.null(treatment)) {
    X_cov <- newX[, setdiff(seq_len(ncol(newX)), c(treat_col, interaction_cols)), drop = FALSE]
    treat_mat <- as.matrix(treatment)
  } else if (length(treat_col) > 0) {
    X_cov <- newX[, setdiff(seq_len(ncol(newX)), c(treat_col, interaction_cols)), drop = FALSE]
    treat_mat <- as.matrix(newX[, treat_col, drop = FALSE])
  } else {
    stop("FindIt requires treatment information in newX or via 'treatment'.", call. = FALSE)
  }

  center <- attr(model, "findit_center")
  scale <- attr(model, "findit_scale")
  Xstd <- .apply_standardization(X_cov, center, scale)

  if (is.null(colnames(Xstd))) {
    colnames(Xstd) <- paste0("Cov", seq_len(ncol(Xstd)))
  }
  if (is.null(colnames(treat_mat))) {
    colnames(treat_mat) <- paste0("treat", seq_len(ncol(treat_mat)))
  }

  start <- model.matrix(~Xstd * treat_mat)
  treat_start <- attr(model, "findit_treat_start")
  treat2 <- start[, treat_start:ncol(start), drop = FALSE]

  pred <- (treat2 %*% model$coefs) / 2 + 0.5
  .clip(as.vector(pred), 0, 1)
}

# =============================================================================
# GLMBoost (mboost)
# =============================================================================

#' Fit GLMBoost model
#' @keywords internal
.fit_glmboost <- function(X, y, treatment, family = "gaussian",
                          weights = NULL, ...) {
  .require_package("mboost", "glmboost")

  X <- .as_matrix(X)
  df <- as.data.frame(X)
  df$.y <- y

  fam <- if (family == "binomial") mboost::Binomial() else mboost::Gaussian()

  model <- mboost::glmboost(
    .y ~ .,
    data = df,
    family = fam,
    weights = weights,
    ...
  )

  attr(model, "family") <- family
  attr(model, "safe_names") <- colnames(df)
  model
}

#' Predict from GLMBoost model
#' @keywords internal
.predict_glmboost <- function(model, newX, treatment = NULL, family = "gaussian") {
  .require_package("mboost", "glmboost")

  newX <- .as_matrix(newX)
  newdf <- as.data.frame(newX)

  preds <- predict(model, newdata = newdf, type = "response")
  as.vector(preds)
}

# =============================================================================
# GLM (baseline)
# =============================================================================

#' Fit GLM model
#' @keywords internal
.fit_glm <- function(X, y, treatment, family = "gaussian", weights = NULL, ...) {
  X <- .as_matrix(X)
  df <- as.data.frame(X)
  df$.y <- y

  fam <- if (family == "binomial") stats::binomial() else stats::gaussian()

  stats::glm(.y ~ .,
             data = df,
             family = fam,
             weights = weights,
             ...)
}

#' Predict from GLM model
#' @keywords internal
.predict_glm <- function(model, newX, treatment = NULL, family = "gaussian") {
  newX <- .as_matrix(newX)
  newdf <- as.data.frame(newX)

  type <- if (family == "binomial") "response" else "response"
  as.vector(stats::predict(model, newdata = newdf, type = type))
}

# =============================================================================
# BART (dbarts)
# =============================================================================

#' Fit BART model
#' @keywords internal
.fit_bart <- function(X, y, treatment, family = "gaussian", weights = NULL, ...) {
  .require_package("dbarts", "bart")

  X <- .as_matrix(X)

  # BART settings
  binary <- family == "binomial"

  # dbarts::bart2 for better interface
  if (binary) {
    model <- dbarts::bart2(
      x.train = X,
      y.train = y,
      n.trees = 200,
      n.samples = 1000,
      n.burn = 250,
      verbose = FALSE,
      ...
    )
  } else {
    model <- dbarts::bart2(
      x.train = X,
      y.train = y,
      n.trees = 200,
      n.samples = 1000,
      n.burn = 250,
      verbose = FALSE,
      ...
    )
  }

  attr(model, "family") <- family
  model
}

#' Predict from BART model
#' @keywords internal
.predict_bart <- function(model, newX, treatment = NULL, family = "gaussian") {
  .require_package("dbarts", "bart")

  newX <- .as_matrix(newX)

  # Get predictions - use mean of posterior samples
  preds <- predict(model, newdata = newX)

  # preds is matrix: n_samples x n_obs for bart2
  # Take column means
  if (is.matrix(preds)) {
    if (nrow(preds) > ncol(preds)) {
      # Samples in rows, obs in cols
      pred_mean <- colMeans(preds)
    } else {
      # Obs in rows, samples in cols
      pred_mean <- rowMeans(preds)
    }
  } else {
    pred_mean <- as.vector(preds)
  }

  # Transform for binary
  if (family == "binomial") {
    pred_mean <- .logistic(pred_mean)
  }

  pred_mean
}

# =============================================================================
# Random Forest (randomForest)
# =============================================================================

#' Fit Random Forest model
#' @keywords internal
.fit_randomforest <- function(X, y, treatment, family = "gaussian",
                              weights = NULL, ...) {
  .require_package("randomForest", "randomForest")

  X <- .as_matrix(X)

  # Convert to data frame with safe names
  df <- as.data.frame(X)
  safe_names <- make.names(colnames(X), unique = TRUE)
  colnames(df) <- safe_names

  if (family == "binomial") {
    # Classification
    y_factor <- factor(y, levels = c(0, 1))

    model <- randomForest::randomForest(
      x = df,
      y = y_factor,
      ntree = 500,
      importance = FALSE,
      ...
    )
  } else {
    # Regression
    model <- randomForest::randomForest(
      x = df,
      y = y,
      ntree = 500,
      importance = FALSE,
      ...
    )
  }

  attr(model, "safe_names") <- safe_names
  attr(model, "family") <- family
  model
}

#' Predict from Random Forest model
#' @keywords internal
.predict_randomforest <- function(model, newX, treatment = NULL, family = "gaussian") {
  .require_package("randomForest", "randomForest")

  newX <- .as_matrix(newX)

  # Use safe names
  safe_names <- attr(model, "safe_names")
  if (!is.null(safe_names) && ncol(newX) == length(safe_names)) {
    colnames(newX) <- safe_names
  }

  newdf <- as.data.frame(newX)

  family <- attr(model, "family") %||% family

  if (family == "binomial") {
    # Get probability of class 1
    preds <- predict(model, newdata = newdf, type = "prob")
    if (is.matrix(preds) && "1" %in% colnames(preds)) {
      as.vector(preds[, "1"])
    } else if (is.matrix(preds)) {
      as.vector(preds[, 2])
    } else {
      as.vector(preds)
    }
  } else {
    as.vector(predict(model, newdata = newdf))
  }
}

# =============================================================================
# KRLS (Kernel Regularized Least Squares)
# =============================================================================

#' Fit KRLS model
#' @keywords internal
.fit_krls <- function(X, y, treatment, family = "gaussian", weights = NULL, ...) {
  .require_package("KRLS", "krls")

  X <- .as_matrix(X)

  # KRLS doesn't support weights directly
  # For binary outcomes, use linear probability model interpretation

  model <- KRLS::krls(
    X = X,
    y = y,
    print.level = 0,
    ...
  )

  attr(model, "family") <- family
  model
}

#' Predict from KRLS model
#' @keywords internal
.predict_krls <- function(model, newX, treatment = NULL, family = "gaussian") {
  .require_package("KRLS", "krls")

  newX <- .as_matrix(newX)

  preds <- predict(model, newdata = newX)

  # Extract fitted values
  if (is.list(preds)) {
    pred_vec <- as.vector(preds$fit)
  } else {
    pred_vec <- as.vector(preds)
  }

  family <- attr(model, "family") %||% family

  if (family == "binomial") {
    # Clip to valid probability range
    pred_vec <- .clip(pred_vec, 0, 1)
  }

  pred_vec
}

# =============================================================================
# SVM (e1071)
# =============================================================================

#' Fit SVM model
#' @keywords internal
.fit_svm <- function(X, y, treatment, family = "gaussian", weights = NULL, ...) {
  .require_package("e1071", "svm")

  X <- .as_matrix(X)

  # Create data frame
  df <- as.data.frame(X)
  safe_names <- make.names(colnames(X), unique = TRUE)
  colnames(df) <- safe_names
  df$.y <- y

  if (family == "binomial") {
    # Classification
    df$.y <- factor(df$.y, levels = c(0, 1))

    model <- e1071::svm(
      .y ~ .,
      data = df,
      type = "C-classification",
      kernel = "radial",
      probability = TRUE,
      ...
    )
  } else {
    # Regression
    model <- e1071::svm(
      .y ~ .,
      data = df,
      type = "eps-regression",
      kernel = "radial",
      ...
    )
  }

  attr(model, "safe_names") <- safe_names
  attr(model, "family") <- family
  model
}

#' Predict from SVM model
#' @keywords internal
.predict_svm <- function(model, newX, treatment = NULL, family = "gaussian") {
  .require_package("e1071", "svm")

  newX <- .as_matrix(newX)

  safe_names <- attr(model, "safe_names")
  if (!is.null(safe_names) && ncol(newX) == length(safe_names)) {
    colnames(newX) <- safe_names
  }

  newdf <- as.data.frame(newX)

  family <- attr(model, "family") %||% family

  if (family == "binomial") {
    preds <- predict(model, newdata = newdf, probability = TRUE)
    probs <- attr(preds, "probabilities")

    if (!is.null(probs) && "1" %in% colnames(probs)) {
      as.vector(probs[, "1"])
    } else if (!is.null(probs)) {
      as.vector(probs[, 2])
    } else {
      # Fallback: convert factor to numeric
      as.numeric(as.character(preds))
    }
  } else {
    as.vector(predict(model, newdata = newdf))
  }
}
