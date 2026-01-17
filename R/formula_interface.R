#' @title Formula Interface for HetEffects
#' @description Parse formulas with pipe syntax for treatment specification.
#' @name formula_interface
#' @keywords internal
NULL

#' Parse formula with pipe syntax
#'
#' Parses formulas of the form \code{outcome ~ covariates | treatment}
#' where the pipe separates covariates from the treatment variable.
#'
#' @param formula Formula object
#' @param data Data frame
#' @return List with outcome, treatment, and covariate names
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' .parse_formula(y ~ x1 + x2 | treat, data)
#' }
.parse_formula <- function(formula, data) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object.", call. = FALSE)
  }

  # Convert to character to find pipe
  fchar <- deparse(formula, width.cutoff = 500)
  fchar <- paste(fchar, collapse = " ")

  # Check for pipe syntax
  has_pipe <- grepl("\\|", fchar)

  if (has_pipe) {
    # Split on pipe
    parts <- strsplit(fchar, "\\|")[[1]]
    if (length(parts) != 2) {
      stop("Formula should have exactly one '|' separating covariates from treatment.",
           call. = FALSE)
    }

    # Parse outcome ~ covariates part
    outcome_covar_part <- trimws(parts[1])
    treatment_part <- trimws(parts[2])

    # Extract outcome and covariates
    outcome_covar_formula <- as.formula(outcome_covar_part)
    parsed <- .parse_simple_formula(outcome_covar_formula, data)

    # Treatment variable
    treatment <- treatment_part
    if (!(treatment %in% names(data))) {
      stop(sprintf("Treatment variable '%s' not found in data.", treatment),
           call. = FALSE)
    }

    list(
      outcome = parsed$outcome,
      treatment = treatment,
      covariates = parsed$covariates,
      formula = outcome_covar_formula
    )

  } else {
    # No pipe - need explicit treatment argument
    parsed <- .parse_simple_formula(formula, data)
    list(
      outcome = parsed$outcome,
      treatment = NULL,  # Must be provided separately
      covariates = parsed$covariates,
      formula = formula
    )
  }
}

#' Parse a simple formula (outcome ~ covariates)
#'
#' @param formula Formula object
#' @param data Data frame
#' @return List with outcome and covariate names
#' @keywords internal
.parse_simple_formula <- function(formula, data) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object.", call. = FALSE)
  }

  # Get terms
  tt <- terms(formula, data = data)

  # Get outcome name (LHS)
  response_var <- all.vars(formula)[1]
  if (attr(tt, "response") == 0) {
    stop("Formula must have a response variable on the left-hand side.", call. = FALSE)
  }

  # Get covariate names from RHS
  # Use model.frame to expand any interactions, etc.
  mf <- model.frame(formula, data = data, na.action = stats::na.pass)
  response_name <- names(mf)[1]

  # Get term labels (covariate names)
  covar_names <- attr(tt, "term.labels")

  # Remove any that are the response
  covar_names <- setdiff(covar_names, response_name)

  # Check they exist in data
  # Handle interaction terms and transformations
  all_vars <- all.vars(formula)
  outcome <- all_vars[1]
  covar_vars <- all_vars[-1]

  missing_vars <- setdiff(covar_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(sprintf("Variables not found in data: %s",
                 paste(missing_vars, collapse = ", ")), call. = FALSE)
  }

  list(
    outcome = outcome,
    covariates = covar_vars,
    terms = tt
  )
}

#' Build model matrix from formula and data
#'
#' @param formula Formula object (without treatment)
#' @param data Data frame
#' @param treatment Treatment vector
#' @param include_interactions Include covariate-treatment interactions
#' @param sparse Return sparse matrix
#' @return Model matrix
#' @keywords internal
.build_model_matrix <- function(formula, data, treatment = NULL,
                                include_interactions = TRUE,
                                sparse = FALSE) {

  # Create model frame and matrix
  mf <- model.frame(formula, data = data, na.action = stats::na.pass)
  y <- model.response(mf)

  # Create model matrix (without intercept by default for penalized methods)
  tt <- terms(formula, data = data)
  X <- model.matrix(tt, data = mf)

  # Remove intercept if present
  intercept_col <- which(colnames(X) == "(Intercept)")
  if (length(intercept_col) > 0) {
    X <- X[, -intercept_col, drop = FALSE]
  }

  # Add treatment interactions if requested
  if (!is.null(treatment) && include_interactions) {
    X <- .make_interaction_matrix(X, treatment, include_main = TRUE)
  } else if (!is.null(treatment)) {
    # Just add treatment as column
    X <- cbind(X, treat = as.numeric(treatment))
  }

  if (sparse && .is_installed("Matrix")) {
    X <- Matrix::Matrix(X, sparse = TRUE)
  }

  list(
    X = X,
    y = y
  )
}

#' Resolve formula vs explicit arguments
#'
#' @param formula Formula (optional)
#' @param data Data frame
#' @param outcome Outcome name (optional)
#' @param treatment Treatment name (optional)
#' @param covariates Covariate names (optional)
#' @return List with resolved outcome, treatment, covariates
#' @keywords internal
.resolve_formula_args <- function(formula = NULL, data, outcome = NULL,
                                  treatment = NULL, covariates = NULL) {

  if (!is.null(formula)) {
    # Use formula interface
    parsed <- .parse_formula(formula, data)

    # Formula-provided values take precedence, but explicit args can override treatment
    resolved_outcome <- parsed$outcome
    resolved_covariates <- parsed$covariates

    # Treatment: formula pipe syntax or explicit arg
    if (!is.null(parsed$treatment)) {
      resolved_treatment <- parsed$treatment
    } else if (!is.null(treatment)) {
      resolved_treatment <- treatment
    } else {
      stop("Treatment must be specified either in formula (y ~ x | treat) or via 'treatment' argument.",
           call. = FALSE)
    }

  } else {
    # Use explicit arguments
    if (is.null(outcome)) {
      stop("Must provide either 'formula' or 'outcome' argument.", call. = FALSE)
    }
    if (is.null(treatment)) {
      stop("Must provide either 'formula' with pipe syntax or 'treatment' argument.",
           call. = FALSE)
    }
    if (is.null(covariates)) {
      stop("Must provide either 'formula' or 'covariates' argument.", call. = FALSE)
    }

    resolved_outcome <- outcome
    resolved_treatment <- treatment
    resolved_covariates <- covariates
  }

  # Validate resolved values
  .validate_data(data, resolved_outcome, resolved_treatment, resolved_covariates)

  list(
    outcome = resolved_outcome,
    treatment = resolved_treatment,
    covariates = resolved_covariates
  )
}
