#' @title Learner Registry and Configuration
#' @description Management of available learners and custom learner registration.
#' @name learner_registry
NULL

# Package-level environment for learner registry
.learner_env <- new.env(parent = emptyenv())
.learner_env$registry <- list()
.learner_env$initialized <- FALSE

#' Initialize the default learner registry
#' @keywords internal
.init_learner_registry <- function() {
  if (.learner_env$initialized) {
    return(invisible())
  }

  # Core learners (glmnet - always available)
  .learner_env$registry$lasso <- list(
    name = "lasso",
    package = "glmnet",
    fit_fn = .fit_lasso,
    predict_fn = .predict_lasso,
    description = "LASSO regression (L1 penalty)",
    category = "core"
  )

  .learner_env$registry$elastic_net_0.75 <- list(
    name = "elastic_net_0.75",
    package = "glmnet",
    fit_fn = function(...) .fit_elastic_net(..., alpha = 0.75),
    predict_fn = .predict_elastic_net,
    description = "Elastic net (alpha = 0.75)",
    category = "core"
  )

  .learner_env$registry$elastic_net_0.5 <- list(
    name = "elastic_net_0.5",
    package = "glmnet",
    fit_fn = function(...) .fit_elastic_net(..., alpha = 0.5),
    predict_fn = .predict_elastic_net,
    description = "Elastic net (alpha = 0.5)",
    category = "core"
  )

  .learner_env$registry$elastic_net_0.25 <- list(
    name = "elastic_net_0.25",
    package = "glmnet",
    fit_fn = function(...) .fit_elastic_net(..., alpha = 0.25),
    predict_fn = .predict_elastic_net,
    description = "Elastic net (alpha = 0.25)",
    category = "core"
  )

  # Original paper learners
  .learner_env$registry$bayesglm <- list(
    name = "bayesglm",
    package = "arm",
    fit_fn = .fit_bayesglm,
    predict_fn = .predict_bayesglm,
    description = "Bayesian GLM with weakly informative priors",
    category = "original"
  )

  .learner_env$registry$bart <- list(
    name = "bart",
    package = "dbarts",
    fit_fn = .fit_bart,
    predict_fn = .predict_bart,
    description = "Bayesian Additive Regression Trees",
    category = "original"
  )

  .learner_env$registry$randomforest <- list(
    name = "randomforest",
    package = "randomForest",
    fit_fn = .fit_randomforest,
    predict_fn = .predict_randomforest,
    description = "Random Forest",
    category = "original"
  )

  .learner_env$registry$krls <- list(
    name = "krls",
    package = "KRLS",
    fit_fn = .fit_krls,
    predict_fn = .predict_krls,
    description = "Kernel Regularized Least Squares",
    category = "original"
  )

  .learner_env$registry$svm <- list(
    name = "svm",
    package = "e1071",
    fit_fn = .fit_svm,
    predict_fn = .predict_svm,
    description = "Support Vector Machine",
    category = "original"
  )

  .learner_env$registry$svm_smo <- list(
    name = "svm_smo",
    package = "e1071",
    fit_fn = .fit_svm,
    predict_fn = .predict_svm,
    description = "SVM (SMO-style, e1071 backend)",
    category = "original"
  )

  .learner_env$registry$findit <- list(
    name = "findit",
    package = "FindIt",
    fit_fn = .fit_findit,
    predict_fn = .predict_findit,
    description = "FindIt (separate penalties for covariates/treatments)",
    category = "original"
  )

  .learner_env$registry$glmboost <- list(
    name = "glmboost",
    package = "mboost",
    fit_fn = .fit_glmboost,
    predict_fn = .predict_glmboost,
    description = "GLMBoost (gradient boosting for GLMs)",
    category = "original"
  )

  .learner_env$registry$glm <- list(
    name = "glm",
    package = "stats",
    fit_fn = .fit_glm,
    predict_fn = .predict_glm,
    description = "Generalized linear model (baseline)",
    category = "original"
  )

  # Modern HTE learners
  .learner_env$registry$causal_forest <- list(
    name = "causal_forest",
    package = "grf",
    fit_fn = .fit_causal_forest,
    predict_fn = .predict_causal_forest,
    description = "Causal Forest (Athey & Wager)",
    category = "modern"
  )

  .learner_env$registry$xgboost <- list(
    name = "xgboost",
    package = "xgboost",
    fit_fn = .fit_xgboost,
    predict_fn = .predict_xgboost,
    description = "Gradient Boosted Trees",
    category = "modern"
  )

  .learner_env$initialized <- TRUE
  invisible()
}

#' Get Available Learners
#'
#' Returns information about learners available in the package.
#'
#' @param check_installed Logical. If TRUE (default), only returns learners
#'   whose required packages are installed.
#' @param category Character. Filter by category: "core", "original", "modern",
#'   or "all" (default).
#'
#' @return A data frame with columns: name, package, description, category,
#'   installed (if check_installed = TRUE).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # See all learners
#' available_learners(check_installed = FALSE)
#'
#' # See only installed learners
#' available_learners()
#'
#' # See only modern HTE learners
#' available_learners(category = "modern")
#' }
available_learners <- function(check_installed = TRUE, category = "all") {
  .init_learner_registry()

  registry <- .learner_env$registry

  df <- data.frame(
    name = sapply(registry, `[[`, "name"),
    package = sapply(registry, `[[`, "package"),
    description = sapply(registry, `[[`, "description"),
    category = sapply(registry, `[[`, "category"),
    stringsAsFactors = FALSE
  )

  if (check_installed) {
    df$installed <- sapply(df$package, .is_installed)
    df <- df[df$installed, ]
    df$installed <- NULL
  }

  if (category != "all") {
    df <- df[df$category == category, ]
  }

  rownames(df) <- NULL
  df
}

#' Register a Custom Learner
#'
#' Adds a custom learner to the registry for use with \code{het_ensemble}.
#'
#' @param name Character. Unique name for the learner.
#' @param fit_fn Function. Fitting function with signature
#'   \code{function(X, y, treatment, family, weights, ...)}.
#' @param predict_fn Function. Prediction function with signature
#'   \code{function(model, newX, treatment, family)}.
#' @param package Character. Name of required package (optional, for documentation).
#' @param description Character. Description of the learner.
#' @param overwrite Logical. If TRUE, overwrite existing learner with same name.
#'
#' @return Invisible NULL.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Register a custom ridge regression learner
#' register_learner(
#'   name = "ridge",
#'   fit_fn = function(X, y, treatment, family, weights, ...) {
#'     glmnet::cv.glmnet(X, y, alpha = 0, family = family, weights = weights)
#'   },
#'   predict_fn = function(model, newX, treatment, family) {
#'     as.vector(predict(model, newX, s = "lambda.min", type = "response"))
#'   },
#'   package = "glmnet",
#'   description = "Ridge regression (L2 penalty)"
#' )
#' }
register_learner <- function(name, fit_fn, predict_fn, package = NULL,
                             description = NULL, overwrite = FALSE) {
  .init_learner_registry()

  if (!is.character(name) || length(name) != 1) {
    stop("'name' must be a single character string.", call. = FALSE)
  }

  if (!is.function(fit_fn)) {
    stop("'fit_fn' must be a function.", call. = FALSE)
  }

  if (!is.function(predict_fn)) {
    stop("'predict_fn' must be a function.", call. = FALSE)
  }

  if (name %in% names(.learner_env$registry) && !overwrite) {
    stop(sprintf("Learner '%s' already exists. Use overwrite = TRUE to replace.",
                 name), call. = FALSE)
  }

  .learner_env$registry[[name]] <- list(
    name = name,
    package = package %||% "custom",
    fit_fn = fit_fn,
    predict_fn = predict_fn,
    description = description %||% "Custom learner",
    category = "custom"
  )

  invisible(NULL)
}

#' Get learner specification by name
#' @param name Learner name
#' @return Learner specification list
#' @keywords internal
.get_learner <- function(name) {
  .init_learner_registry()

  if (!(name %in% names(.learner_env$registry))) {
    stop(sprintf("Unknown learner: '%s'. See available_learners().", name),
         call. = FALSE)
  }

  .learner_env$registry[[name]]
}

#' Get learner specifications for a vector of names
#' @param names Character vector of learner names
#' @param check_installed Check if packages are installed
#' @return Named list of learner specifications
#' @keywords internal
.get_learners <- function(names, check_installed = TRUE) {
  .init_learner_registry()

  learners <- lapply(names, function(nm) {
    spec <- .get_learner(nm)

    if (check_installed && !.is_installed(spec$package)) {
      warning(sprintf("Package '%s' required for learner '%s' is not installed. Skipping.",
                      spec$package, nm), call. = FALSE)
      return(NULL)
    }

    spec
  })

  names(learners) <- names
  learners <- learners[!sapply(learners, is.null)]
  learners
}

#' Auto-detect available learners
#'
#' Returns learner names for all learners with installed packages.
#'
#' @return Character vector of learner names
#' @keywords internal
.auto_detect_learners <- function() {
  .init_learner_registry()

  registry <- .learner_env$registry

  installed <- sapply(registry, function(spec) {
    .is_installed(spec$package)
  })

  names(registry)[installed]
}

#' Check learner packages and warn about missing ones
#' @param learners Character vector of learner names
#' @return Character vector of available learner names
#' @keywords internal
.check_learner_packages <- function(learners) {
  .init_learner_registry()

  available <- character(0)
  missing <- character(0)

  for (nm in learners) {
    spec <- .get_learner(nm)
    if (.is_installed(spec$package)) {
      available <- c(available, nm)
    } else {
      missing <- c(missing, sprintf("%s (requires %s)", nm, spec$package))
    }
  }

  if (length(missing) > 0) {
    warning(sprintf("Some learners unavailable due to missing packages: %s",
                    paste(missing, collapse = ", ")), call. = FALSE)
  }

  available
}

#' Null coalescing operator
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
