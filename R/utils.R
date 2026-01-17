#' @title Utility Functions for HetEffects
#' @description Internal helper functions used throughout the package.
#' @name utils
#' @keywords internal
NULL

#' Logistic (sigmoid) function
#' @param x Numeric vector
#' @return Numeric vector of probabilities
#' @keywords internal
.logistic <- function(x) {
  1 / (1 + exp(-x))
}

#' Inverse logistic (logit) function
#' @param p Numeric vector of probabilities
#' @return Numeric vector
#' @keywords internal
.logit <- function(p) {
  log(p / (1 - p))
}

#' Standardize a numeric vector
#' @param x Numeric vector
#' @param na.rm Logical, remove NA values
#' @return Standardized vector (mean 0, sd 1)
#' @keywords internal
.standardize <- function(x, na.rm = TRUE) {
  m <- mean(x, na.rm = na.rm)
  s <- sd(x, na.rm = na.rm)
  if (s == 0 || is.na(s)) {
    return(x - m)
  }

  (x - m) / s
}

#' Standardize columns of a matrix
#' @param X Numeric matrix
#' @param center Logical, center columns
#' @param scale Logical, scale columns
#' @return List with standardized matrix and transformation parameters
#' @keywords internal
.standardize_matrix <- function(X, center = TRUE, scale = TRUE) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  means <- colMeans(X, na.rm = TRUE)
  sds <- apply(X, 2, sd, na.rm = TRUE)
  sds[sds == 0 | is.na(sds)] <- 1

  if (center) {
    X <- sweep(X, 2, means, "-")
  } else {
    means <- rep(0, p)
  }

  if (scale) {
    X <- sweep(X, 2, sds, "/")
  } else {
    sds <- rep(1, p)
  }

  list(
    X = X,
    center = means,
    scale = sds
  )
}

#' Apply standardization parameters to new data
#' @param X Numeric matrix
#' @param center Vector of means
#' @param scale Vector of sds
#' @return Standardized matrix
#' @keywords internal
.apply_standardization <- function(X, center, scale) {
  X <- as.matrix(X)
  X <- sweep(X, 2, center, "-")
  X <- sweep(X, 2, scale, "/")
  X
}

#' Create interaction matrix between covariates and treatment
#' @param X Covariate matrix
#' @param treatment Treatment vector
#' @param include_main Logical, include main effects
#' @return Interaction matrix
#' @keywords internal
.make_interaction_matrix <- function(X, treatment, include_main = TRUE) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  # Convert treatment to numeric if factor

  if (is.data.frame(treatment)) {
    treatment <- as.matrix(treatment)
  }

  if (is.matrix(treatment)) {
    treat_mat <- treatment
  } else {
    if (is.factor(treatment)) {
      treatment <- as.numeric(treatment) - 1
    }
    treat_mat <- as.numeric(treatment)
  }

  # Create interaction terms
  if (is.matrix(treat_mat)) {
    X_interact <- do.call(cbind, lapply(seq_len(ncol(treat_mat)), function(j) {
      X * treat_mat[, j]
    }))
  } else {
    X_interact <- X * treat_mat
  }

  # Set column names
  orig_names <- colnames(X)
  if (is.null(orig_names)) {
    orig_names <- paste0("X", seq_len(p))
  }
  if (is.matrix(treat_mat)) {
    treat_names <- colnames(treat_mat)
    if (is.null(treat_names)) {
      treat_names <- paste0("treat", seq_len(ncol(treat_mat)))
    }
    colnames(X_interact) <- unlist(lapply(treat_names, function(tn) {
      paste0(orig_names, ":", tn)
    }))
  } else {
    colnames(X_interact) <- paste0(orig_names, ":treat")
  }

  if (include_main) {
    if (is.matrix(treat_mat)) {
      if (is.null(colnames(treat_mat))) {
        colnames(treat_mat) <- paste0("treat", seq_len(ncol(treat_mat)))
      }
      result <- cbind(X, treat_mat, X_interact)
    } else {
      result <- cbind(X, treat = treat_mat, X_interact)
    }
  } else {
    if (is.matrix(treat_mat)) {
      result <- cbind(treat_mat, X_interact)
    } else {
      result <- cbind(treat = treat_mat, X_interact)
    }
  }

  result
}

#' Safe matrix conversion
#' @param x Object to convert
#' @return Matrix
#' @keywords internal
.as_matrix <- function(x) {
  if (inherits(x, "Matrix")) {
    as.matrix(x)
  } else if (is.data.frame(x)) {
    as.matrix(x)
  } else {
    as.matrix(x)
  }
}

#' Check if a package is installed
#' @param pkg Package name
#' @return Logical
#' @keywords internal
.is_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

#' Safely load a package
#' @param pkg Package name
#' @param fn Function name to check
#' @return Logical indicating success
#' @keywords internal
.require_package <- function(pkg, fn = NULL) {
  if (!.is_installed(pkg)) {
    if (!is.null(fn)) {
      stop(sprintf("Package '%s' is required for function '%s'. Please install it.", pkg, fn),
           call. = FALSE)
    } else {
      stop(sprintf("Package '%s' is required. Please install it.", pkg), call. = FALSE)
    }
  }
  TRUE
}

#' Weighted mean with NA handling
#' @param x Numeric vector
#' @param w Weights
#' @param na.rm Remove NAs
#' @return Weighted mean
#' @keywords internal
.weighted_mean <- function(x, w = NULL, na.rm = TRUE) {
  if (is.null(w)) {
    return(mean(x, na.rm = na.rm))
  }

  if (na.rm) {
    idx <- !is.na(x) & !is.na(w)
    x <- x[idx]
    w <- w[idx]
  }

  if (length(x) == 0) return(NA_real_)
  sum(x * w) / sum(w)
}

#' Weighted variance
#' @param x Numeric vector
#' @param w Weights
#' @param na.rm Remove NAs
#' @return Weighted variance
#' @keywords internal
.weighted_var <- function(x, w = NULL, na.rm = TRUE) {
  if (is.null(w)) {
    return(var(x, na.rm = na.rm))
  }

  if (na.rm) {
    idx <- !is.na(x) & !is.na(w)
    x <- x[idx]
    w <- w[idx]
  }


  if (length(x) <= 1) return(NA_real_)

  # Reliability weights (frequency weights)
  wm <- .weighted_mean(x, w, na.rm = FALSE)
  n <- length(x)
  sum(w * (x - wm)^2) / (sum(w) - sum(w)^2 / sum(w^2))
}

#' Clip values to a range
#' @param x Numeric vector
#' @param lower Lower bound
#' @param upper Upper bound
#' @return Clipped vector
#' @keywords internal
.clip <- function(x, lower = -Inf, upper = Inf) {
  pmax(pmin(x, upper), lower)
}

#' Check that treatment is binary
#' @param treatment Treatment vector
#' @return Logical
#' @keywords internal
.check_binary_treatment <- function(treatment) {
  if (is.data.frame(treatment)) {
    treatment <- as.matrix(treatment)
  }
  if (is.matrix(treatment)) {
    apply(treatment, 2, function(col) {
      vals <- unique(na.omit(col))
      if (length(vals) != 2) {
        stop("Each treatment column must be binary (exactly 2 unique values).",
             call. = FALSE)
      }
      TRUE
    })
    return(TRUE)
  }

  vals <- unique(na.omit(treatment))
  if (length(vals) != 2) {
    stop("Treatment must be binary (exactly 2 unique values).", call. = FALSE)
  }
  TRUE
}

#' Convert treatment to 0/1 coding
#' @param treatment Treatment vector
#' @param control Control value
#' @param treated Treated value
#' @return Numeric 0/1 vector
#' @keywords internal
.recode_treatment <- function(treatment, control = NULL, treated = NULL) {
  if (is.data.frame(treatment)) {
    treatment <- as.matrix(treatment)
  }
  if (is.matrix(treatment)) {
    recoded <- apply(treatment, 2, function(col) {
      vals <- sort(unique(na.omit(col)))
      if (length(vals) != 2) {
        stop("Each treatment column must be binary (exactly 2 unique values).",
             call. = FALSE)
      }
      result <- rep(NA_real_, length(col))
      result[col == vals[1]] <- 0
      result[col == vals[2]] <- 1
      result
    })
    colnames(recoded) <- colnames(treatment)
    return(recoded)
  }

  vals <- sort(unique(na.omit(treatment)))

  if (is.null(control) && is.null(treated)) {
    # Auto-detect: smaller value is control
    control <- vals[1]
    treated <- vals[2]
  } else if (is.null(control)) {
    control <- setdiff(vals, treated)
  } else if (is.null(treated)) {
    treated <- setdiff(vals, control)
  }

  result <- rep(NA_real_, length(treatment))
  result[treatment == control] <- 0
  result[treatment == treated] <- 1
  result
}

#' Create a progress message
#' @param msg Message
#' @param verbose Logical
#' @keywords internal
.message <- function(msg, verbose = TRUE) {
  if (verbose) {
    message(msg)
  }
}

#' Safe try-catch wrapper
#' @param expr Expression to evaluate
#' @param on_error How to handle errors
#' @param learner_name Name of learner (for messages)
#' @return Result or NULL on error
#' @keywords internal
.safe_try <- function(expr, on_error = c("warn", "stop", "ignore"),
                      learner_name = "learner") {
  on_error <- match.arg(on_error)

  result <- tryCatch(
    expr,
    error = function(e) {
      if (on_error == "stop") {
        stop(sprintf("Error in %s: %s", learner_name, e$message), call. = FALSE)
      } else if (on_error == "warn") {
        warning(sprintf("Error in %s: %s. Skipping.", learner_name, e$message),
                call. = FALSE)
      }
      NULL
    }
  )

  result
}

#' Validate input data
#' @param data Data frame
#' @param outcome Outcome column
#' @param treatment Treatment column
#' @param covariates Covariate columns
#' @return Validated data components
#' @keywords internal
.validate_data <- function(data, outcome, treatment, covariates) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }

  # Check outcome

  if (!(outcome %in% names(data))) {
    stop(sprintf("Outcome variable '%s' not found in data.", outcome), call. = FALSE)
  }
  y <- data[[outcome]]

  # Check treatment
  if (length(treatment) > 1) {
    missing_treat <- setdiff(treatment, names(data))
    if (length(missing_treat) > 0) {
      stop(sprintf("Treatment variables not found in data: %s",
                   paste(missing_treat, collapse = ", ")), call. = FALSE)
    }
    treat <- data[, treatment, drop = FALSE]
  } else {
    if (!(treatment %in% names(data))) {
      stop(sprintf("Treatment variable '%s' not found in data.", treatment), call. = FALSE)
    }
    treat <- data[[treatment]]
  }
  .check_binary_treatment(treat)

  # Check covariates
  missing_covars <- setdiff(covariates, names(data))
  if (length(missing_covars) > 0) {
    stop(sprintf("Covariates not found in data: %s",
                 paste(missing_covars, collapse = ", ")), call. = FALSE)
  }
  X <- data[, covariates, drop = FALSE]

  list(
    y = y,
    treatment = treat,
    X = X,
    n = nrow(data),
    p = length(covariates)
  )
}

#' Get default family based on outcome
#' @param y Outcome vector
#' @return Character, "binomial" or "gaussian"
#' @keywords internal
.detect_family <- function(y) {
  vals <- unique(na.omit(y))
  if (length(vals) == 2 && all(vals %in% c(0, 1))) {
    "binomial"
  } else {
    "gaussian"
  }
}
