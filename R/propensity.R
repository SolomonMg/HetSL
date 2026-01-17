#' @title Propensity Score Estimation
#' @description Estimate propensity scores for observational studies.
#' @name propensity
#' @keywords internal
NULL

#' Estimate Propensity Scores
#'
#' @param formula Formula of the form \code{treatment ~ covariates}.
#' @param data Data frame.
#' @param method Method for estimation: "logit", "gbm", or "rf".
#' @param weights Optional observation weights.
#' @param ... Additional arguments passed to the underlying learner.
#'
#' @return A list with elements:
#' \describe{
#'   \item{propensity}{Estimated propensity scores.}
#'   \item{model}{Fitted model object.}
#'   \item{method}{Method used.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ps <- estimate_propensity(treat ~ x1 + x2, data = mydata)
#' }
estimate_propensity <- function(formula, data, method = c("logit", "gbm", "rf"),
                                weights = NULL, ...) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }

  method <- match.arg(method)

  # Validate treatment is binary
  treat_name <- all.vars(formula)[1]
  if (!(treat_name %in% names(data))) {
    stop(sprintf("Treatment variable '%s' not found in data.", treat_name), call. = FALSE)
  }
  .check_binary_treatment(data[[treat_name]])

  if (method == "logit") {
    if (is.null(weights)) {
      model <- stats::glm(
        formula = formula,
        data = data,
        family = stats::binomial(),
        ...
      )
    } else {
      model <- stats::glm(
        formula = formula,
        data = data,
        family = stats::binomial(),
        weights = weights,
        ...
      )
    }
    propensity <- stats::predict(model, type = "response")

  } else if (method == "gbm") {
    if (!.is_installed("gbm")) {
      stop("Package 'gbm' is required for method = 'gbm'. Please install it.",
           call. = FALSE)
    }

    mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
    y <- stats::model.response(mf)
    X <- stats::model.matrix(formula, data = mf)

    if (!is.null(weights)) {
      w <- weights
    } else {
      w <- rep(1, length(y))
    }

    model <- gbm::gbm.fit(
      x = X,
      y = y,
      distribution = "bernoulli",
      weights = w,
      verbose = FALSE,
      ...
    )

    propensity <- gbm::predict.gbm(model, newdata = X, type = "response")

  } else {
    if (!.is_installed("randomForest")) {
      stop("Package 'randomForest' is required for method = 'rf'. Please install it.",
           call. = FALSE)
    }

    mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
    y <- stats::model.response(mf)
    X <- stats::model.matrix(formula, data = mf)

    y_factor <- factor(y, levels = c(0, 1))

    model <- randomForest::randomForest(
      x = X,
      y = y_factor,
      ntree = 500,
      importance = FALSE,
      ...
    )

    prob <- stats::predict(model, newdata = X, type = "prob")
    propensity <- if (is.matrix(prob) && "1" %in% colnames(prob)) {
      prob[, "1"]
    } else {
      prob[, 2]
    }
  }

  list(
    propensity = as.vector(propensity),
    model = model,
    method = method
  )
}
