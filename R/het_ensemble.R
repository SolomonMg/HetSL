#' @title Fit Heterogeneous Treatment Effects Ensemble
#' @description Main function for fitting ensemble models for heterogeneous treatment effects.
#' @name het_ensemble
NULL

#' Fit Heterogeneous Treatment Effects Ensemble
#'
#' Fits an ensemble of machine learning models to estimate heterogeneous
#' treatment effects, following Grimmer, Messing & Westwood (2017).
#'
#' @param formula Formula specifying the model. Can use pipe syntax:
#'   \code{outcome ~ covariates | treatment} where covariates and treatment
#'   are separated by \code{|}. Alternatively, specify outcome, treatment,
#'   and covariates explicitly.
#' @param data A data frame containing the variables in the model.
#' @param outcome Character. Name of the outcome variable (alternative to formula).
#' @param treatment Character. Name of the treatment variable (alternative to formula).
#' @param covariates Character vector. Names of covariate variables (alternative to formula).
#' @param learners Character vector of learner names to use, or NULL to auto-detect
#'   available learners based on installed packages. See \code{\link{available_learners}}.
#' @param nfolds Integer. Number of cross-validation folds (default 10).
#' @param family Character. Either "binomial" for binary outcomes or "gaussian"
#'   for continuous outcomes. If NULL, auto-detected from outcome.
#' @param weights Numeric vector of observation weights (e.g., sampling weights).
#' @param propensity Either a character string naming a propensity score column
#'   in data, or a formula for propensity score estimation. Used for
#'   observational studies.
#' @param parallel Logical. Use parallel processing via future/furrr.
#' @param on_error Character. How to handle learner failures: "warn" (default),
#'   "stop", or "ignore".
#' @param seed Integer. Random seed for reproducibility.
#' @param verbose Logical. Print progress messages.
#'
#' @return An object of class \code{het_ensemble} with components:
#' \describe{
#'   \item{weights}{Named vector of ensemble weights}
#'   \item{models}{List of fitted model objects}
#'   \item{cv_predictions}{Matrix of cross-validation predictions}
#'   \item{cv_performance}{Data frame of CV performance metrics}
#'   \item{call}{The matched call}
#'   \item{data}{List with y, X, treatment vectors}
#'   \item{formula}{Formula used (if any)}
#'   \item{family}{Outcome family}
#'   \item{learner_names}{Names of learners used}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Using formula interface with pipe syntax
#' fit <- het_ensemble(
#'   outcome ~ age + female + education | treatment,
#'   data = mydata,
#'   nfolds = 10
#' )
#'
#' # Using explicit arguments
#' fit <- het_ensemble(
#'   data = mydata,
#'   outcome = "y",
#'   treatment = "treat",
#'   covariates = c("age", "female", "education")
#' )
#'
#' # With specific learners
#' fit <- het_ensemble(
#'   outcome ~ x1 + x2 | treat,
#'   data = mydata,
#'   learners = c("lasso", "randomforest", "causal_forest")
#' )
#'
#' # Observational study with propensity weighting
#' fit <- het_ensemble(
#'   outcome ~ x1 + x2 | treat,
#'   data = mydata,
#'   propensity = "ps_column"  # Pre-estimated propensity scores
#' )
#' }
het_ensemble <- function(formula = NULL,
                         data,
                         outcome = NULL,
                         treatment = NULL,
                         covariates = NULL,
                         learners = NULL,
                         nfolds = 10L,
                         family = NULL,
                         weights = NULL,
                         propensity = NULL,
                         parallel = FALSE,
                         on_error = c("warn", "stop", "ignore"),
                         seed = NULL,
                         verbose = TRUE) {


  # Capture call
  cl <- match.call()
  on_error <- match.arg(on_error)

  # Set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Validate data
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }

  # Resolve formula vs explicit arguments
  resolved <- .resolve_formula_args(
    formula = formula,
    data = data,
    outcome = outcome,
    treatment = treatment,
    covariates = covariates
  )

  outcome_var <- resolved$outcome
  treatment_var <- resolved$treatment
  covariate_vars <- resolved$covariates

  .message(sprintf("Outcome: %s", outcome_var), verbose)
  .message(sprintf("Treatment: %s", treatment_var), verbose)
  .message(sprintf("Covariates: %s", paste(covariate_vars, collapse = ", ")), verbose)

  # Extract data vectors
  y <- data[[outcome_var]]
  treat <- if (length(treatment_var) > 1) {
    data[, treatment_var, drop = FALSE]
  } else {
    data[[treatment_var]]
  }
  X <- data[, covariate_vars, drop = FALSE]

  # Handle missing values
  complete_idx <- complete.cases(y, treat, X)
  if (sum(!complete_idx) > 0) {
    warning(sprintf("Removing %d observations with missing values.",
                    sum(!complete_idx)), call. = FALSE)
    y <- y[complete_idx]
    treat <- treat[complete_idx]
    X <- X[complete_idx, , drop = FALSE]
    if (!is.null(weights)) weights <- weights[complete_idx]
  }

  n <- length(y)
  .message(sprintf("Sample size: %d", n), verbose)

  # Recode treatment to 0/1
  treat_numeric <- .recode_treatment(treat)

  # Detect or validate family
  if (is.null(family)) {
    family <- .detect_family(y)
    .message(sprintf("Auto-detected family: %s", family), verbose)
  }
  family <- match.arg(family, c("binomial", "gaussian"))

  # Handle propensity scores for observational studies
  propensity_weights <- NULL
  if (!is.null(propensity)) {
    if (is.matrix(treat_numeric)) {
      stop("Propensity weighting is only supported for binary treatments.",
           call. = FALSE)
    }
    propensity_weights <- .handle_propensity(propensity, data[complete_idx, ],
                                             treat_numeric, verbose)
    # Combine with observation weights
    if (!is.null(weights)) {
      weights <- weights * propensity_weights
    } else {
      weights <- propensity_weights
    }
  }

  # Build model matrix with interactions
  X_matrix <- .as_matrix(X)
  X_full <- .make_interaction_matrix(X_matrix, treat_numeric, include_main = TRUE)

  .message(sprintf("Design matrix: %d x %d", nrow(X_full), ncol(X_full)), verbose)

  # Determine learners
  if (is.null(learners)) {
    learners <- .auto_detect_learners()
    .message(sprintf("Auto-detected %d available learners.", length(learners)), verbose)
  } else {
    learners <- .check_learner_packages(learners)
  }

  if (length(learners) == 0) {
    stop("No learners available. Install required packages or check learner names.",
         call. = FALSE)
  }

  .message(sprintf("Using learners: %s", paste(learners, collapse = ", ")), verbose)

  # Get learner specifications
  learner_specs <- .get_learners(learners, check_installed = TRUE)

  # Set up parallel processing
  if (parallel) {
    if (.is_installed("future") && .is_installed("furrr")) {
      .message("Using parallel processing via future/furrr.", verbose)
    } else {
      warning("Parallel processing requested but future/furrr not installed. Using sequential.",
              call. = FALSE)
      parallel <- FALSE
    }
  }

  # Run cross-validation
  cv_results <- .run_cv_all_learners(
    X = X_full,
    y = y,
    treatment = treat_numeric,
    learners = learner_specs,
    nfolds = nfolds,
    family = family,
    weights = weights,
    on_error = on_error,
    parallel = parallel,
    seed = seed,
    verbose = verbose
  )

  # Optimize ensemble weights
  .message("Optimizing ensemble weights...", verbose)
  ensemble_weights <- .optimize_weights(
    cv_predictions = cv_results$cv_predictions,
    y = y,
    family = family,
    weights = weights
  )

  .message("Fitting final models...", verbose)

  # Fit final models on full data
  final_models <- lapply(cv_results$successful_learners, function(lname) {
    lspec <- learner_specs[[lname]]

    .safe_try(
      lspec$fit_fn(X_full, y, treat_numeric, family = family, weights = weights),
      on_error = on_error,
      learner_name = lname
    )
  })
  names(final_models) <- cv_results$successful_learners

  # Remove any failed final fits
  failed_final <- sapply(final_models, is.null)
  if (any(failed_final)) {
    warning(sprintf("Final fit failed for: %s",
                    paste(names(final_models)[failed_final], collapse = ", ")),
            call. = FALSE)
    final_models <- final_models[!failed_final]
    ensemble_weights <- ensemble_weights[names(final_models)]
    ensemble_weights <- ensemble_weights / sum(ensemble_weights)
    cv_results$cv_predictions <- cv_results$cv_predictions[, names(final_models), drop = FALSE]
  }

  # Compute CV performance metrics
  cv_performance <- .learner_performance(
    cv_results$cv_predictions, y, family, weights
  )
  cv_performance$ensemble_weight <- ensemble_weights[cv_performance$learner]

  .message("Done.", verbose)

  # Build result object
  result <- list(
    weights = ensemble_weights,
    models = final_models,
    cv_predictions = cv_results$cv_predictions,
    cv_performance = cv_performance,
    call = cl,
    data = list(
      y = y,
      X = X,
      X_full = X_full,
      treatment = treat_numeric,
      treatment_original = treat,
      n = n,
      weights = weights,
      propensity_weights = propensity_weights,
      complete_idx = complete_idx
    ),
    formula = formula,
    family = family,
    learner_names = names(final_models),
    outcome_var = outcome_var,
    treatment_var = treatment_var,
    covariate_vars = covariate_vars,
    nfolds = nfolds,
    folds = cv_results$folds
  )

  class(result) <- "het_ensemble"
  result
}

#' Handle propensity score specification
#' @keywords internal
.handle_propensity <- function(propensity, data, treatment, verbose = TRUE) {
  if (is.character(propensity) && length(propensity) == 1) {
    # Column name
    if (!(propensity %in% names(data))) {
      stop(sprintf("Propensity column '%s' not found in data.", propensity),
           call. = FALSE)
    }
    ps <- data[[propensity]]
    .message("Using pre-estimated propensity scores.", verbose)
  } else if (inherits(propensity, "formula")) {
    # Estimate propensity scores
    .message("Estimating propensity scores...", verbose)
    ps <- estimate_propensity(propensity, data, method = "logit")$propensity
  } else {
    stop("'propensity' must be a column name or formula.", call. = FALSE)
  }

  # Compute IPW weights
  # w = T/p + (1-T)/(1-p) for ATE
  # Stabilized weights: w = T * mean(T)/p + (1-T) * mean(1-T)/(1-p)
  ps <- .clip(ps, 0.01, 0.99)  # Truncate extreme values
  p_treat <- mean(treatment)

  ipw <- treatment / ps + (1 - treatment) / (1 - ps)

  # Normalize to sum to n
  ipw <- ipw * length(ipw) / sum(ipw)

  ipw
}

#' Check if object is het_ensemble
#' @param x Object to check
#' @return Logical
#' @keywords internal
is.het_ensemble <- function(x) {
  inherits(x, "het_ensemble")
}
