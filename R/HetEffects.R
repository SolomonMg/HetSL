# Fits a super learned ensamble to estimate heterogenous treatment effects.
# Key decisions and TODOs:
#   - Do we want users to be able to specify newX themselves?
#     - Do we want to force users to do this themselves,
#       and just provide helper functions?
#   - How to handle continuous variables?
#     - If estimating MCATEs, do we estimate at the mean?
#     - If not, do we fix at the mean?
#   - Still need to process and format results, display for the user.
#     We can base this on previous code.
#   - Test, test, test.
#   -
#   -
###############################################################################

#' Estimating Heterogenous Treatments and Effects using Super Learning
#'
#' The 'HetEffects' command uses the method outlined in Grimmer, Messing
#' and Westwood (2014) to estimate heterogenous treatments and heterogenous
#' treatment effects.  This method uses an ensemble of constituent methods,
#' weighted based on cross-validated model performance.
#' An ensemble of constituent models are fit based on user-supplied formulas,
#' weighted according to cross-validated model performance, then used to
#' predict the expected marginal conditional averages.  Then the treatment
#' and control are differenced, conditional on key user-specified covariates
#' to produce Marginal Conditional Average Treatment Effects (MCATEs).
#'
#'
#' @usage
#' HetEffects(
#'    formula,
#'    data,
#'    treatments,
#'    subset,
#'    weights = NULL,
#'    family = gaussian(),
#'    SL.library,
#'    method = 'method.NNLS',
#'    id = NULL,
#'    verbose = TRUE,
#'    control = list(),
#'    cvControl = list(),
#'    bootstrap = TRUE,
#'    cores = 1,
#'    R = 2,
#'    na.action = NULL, ...)
#'
#' @param formula The formula used to estimate hetergenous effects
#' @param data The data to be used.  Data MUST be in XYZ format.
#'        Data MUST NOT contain NAs.
#' @param treatments Specify which of columns in the data are treatments.
#' @param subset The subset of variables for which to return MCATEs.
#'        Setting to NULL will return all variables in the formula.
#' @param weights observation weights to be used to fit models.
#'        Verify SuperLearner can utilize weights in each model.
#' @param Family link function for models. Currently either gaussian() or binomial().
#' @param SL.library library of algorithms
#' @param method Loss function and model to estimate ensemble weights.
#'        Default is method.NNLS, method.NNloglik, or custom
#'        (see SuperLearner::create.method.template()).
#' @param id cluster id
#' @param verbose View progress of estimation? Defaults to TRUE.
#' @param control Optional controls for the SuperLearner package.
#' @param cvControl List for CV control, see SuperLearner documentation.
#' @param bootstrap Bootstrap estimates?  Defaults to FALSE.
#' @param cores Number of cores to use if using doMC for parallelization.
#' @param R Number of bootstrap replicate samples to estimate if bootstrap == TRUE.
#' @param na.action How to handle NA values.
#' @param \dots additional parameters to be passed to the plot
#'
#' @return if bootstrap is FALSE, returns a list of model weights and
#'         predictions for each model, and the corresponding covariate
#'         combinations, effectX.  If bootstrap is TRUE,
#'         returns a list of the above (of length 1:R) for each boostrap
#'         replicate, and the corresponding covariate
#'         combinations, effectX
#' @export
#' @author Solomon Messing and Sean Westwood and Justin Grimmer


HetEffects <- function(formula,
    data,
    X_intercept = FALSE,
    treatments,
    max_covariates = 1e10,
    subset,
    weights = NULL,
    family = gaussian(),
    SL.library,
    method = 'method.NNLS',
    id = NULL,
    verbose = TRUE,
    control = list(),
    cvControl = list(),
    bootstrap = FALSE,
    cores = 1,
    R = 2,
    na.action = NULL, ...){
  # TODO: allow people to pass covariate matrix + weights?
  # TODO: Implement bootstrapping?
  # TODO: Consider setting tuning params on the fly
    # glmnet::alpha
    # knn::k

  require('SuperLearner')

  # vector for treatment column names. Check to make sure all is well:
  if (!is.character(treatments)) stop("'treatments' must be a character scalar/vector")
  if (!is.vector(treatments)) stop("'treatments' must be a character scalar/vector")
  for (i in 1:length(treatments)){
    if (length(grep(treatments[i], formula, fixed=TRUE))==0) {
      stop("All treatments must be in the formula")
    }
  }
  
  if (X_intercept==FALSE & length(grep("0", formula[3])) == 0){
    stop("Formula must contain '0 +' unless you're including an intercept")
  }
  
  # TODO: implement with Matrix::sparse.model.matrix
  
  # Parse formula, data into model matrix, Y
  X_frame <- model.frame(formula, data=data)
  names(X_frame)
  toInclude <- complete.cases(X_frame)
  
  Xmat <- model.matrix(formula, data=X_frame)
  Xmat <- Xmat[toInclude,]
  Y <- X_frame[[1]]
  Y <- Y[toInclude]
  colnames(Xmat)
  Xmat <- data.frame(Xmat)
  
  # Check to make sure there are no NAs in Xmat values
  anyNA = function(x) any(is.na(x))
  if (any(apply(Xmat, 2, anyNA))){
    stop("NAs remain, remove or recode before execution")
  }
  
  # add character separation between different categories for same variable
  # in design matrix:
  formula_matrix_char_sep <- function(formula, Xmat){
    group_vars <- unique(unlist(strsplit(as.character(formula[3]), 
                                         split=" \\+ | \\* ")))
    for(i in 1:length(group_vars)){
      var_idx <- grep(pattern = group_vars[i], x = colnames(Xmat))
      colnames(Xmat)[var_idx] <- gsub(group_vars[i], 
                  paste0(group_vars[i], "___"), colnames(Xmat)[var_idx])
    }
    Xmat
  }

  Xmat <- formula_matrix_char_sep(formula, Xmat)
  names(Xmat)
  
  # Take all covariates, all treatments, interactions
  # rename treatment columns T_whatever so that we can grep
  # them in FindIt.

  # Find treatment cols
  find_treatment_cols <- function(treatments, Xmat){
    treat_columns = c()
    for(i in 1:length(treatments)){
      tTreat_columns <- c(treat_columns, grep(treatments[i],
                                             colnames(Xmat))) 
    }
    treat_columns = sort(unique(Treat_columns))
    treat_columns
  }
  
  treat_columns <- find_treatment_cols(treatments, Xmat)
  
  T_var_names <- colnames(Xmat)[treat_columns]
  cov_var_names <- colnames(Xmat)[-treat_columns]
  
  # Mark w/ T_
  colnames(Xmat)[treat_columns] <- paste0("T_", colnames((Xmat))[treat_columns])
  # attr(Xmat, "treat_columns") <- treat_columns
  
  # Make sure there are no numeric variables
  col_numeric <- apply(Xmat, 2, function(x) length(unique(x)) > 20 )
  if (any(col_numeric)) stop(
    "Columns in design matrix must be discrete (20 unique values or fewer).")

  # Generate matrix for prediction:
  
  # veclist = apply(Xmat, 2, unique)
  veclist = apply(X_frame[,-1], 2, unique)
  # veclist = data.frame(veclist)
  veclist = lapply(veclist, na.omit)
  # veclist = lapply(veclist, as.vector)
  
  # Make sure there are not too many factor combinations
  n_cov_combinations <- prod(unlist(lapply(veclist, length)))
  if (n_cov_combinations > max_covariates) {
    stop("Too many factor combinations to estimate effects")
  }
  
  # Create matrix for effects:
  # newX = expand.grid(veclist, stringsAsFactors=FALSE)
  newX = expand.grid(veclist)
  
  # for every var in newX, relevel according to X_frame
  same_factor_lvls <- function(old_x, new_x){
    x <- factor(new_x, levels = levels(old_x))
    x
  }
  
  for (col in names(newX)){
    newX[[col]] <- same_factor_lvls(X_frame[[col]], newX[[col]])
    
  }
  
  # Put together model matrix for newX
  formula_rhs <- as.formula(paste0( "~", formula[3]))
  newXmat <- model.matrix(formula_rhs, data=newX)
  newXmat <- formula_matrix_char_sep(formula, newXmat)
  colnames(newXmat)
  
  # Add "T_" to varnames
  treat_columns <- find_treatment_cols(treatments, newXmat)
  colnames(newXmat)[treat_columns] <- paste0("T_", colnames((newXmat))[treat_columns])
  
  newXmatb <- data.frame(newXmat)
  newXmatb <- newXmatb[,colnames(Xmatb)]
  
  # Put everything in mod matrix:
  Xmatb <- data.frame(Xmat)
  
  # Fix colnames (data.frame subsitutes " " with ".")
  # colnames(Xmatb) <- colnames(Xmat)
  setdiff(names(Xmatb), names(newXmatb))
  setdiff(names(newXmatb), names(Xmatb))
  
  attr(Xmatb, "treatments") <- treatments
  
  # No longer needed? 
  # newXb <- data.frame(model.matrix(~ -1 +., newX))
  # names(newXb) == names(Xmatb)
  
  if(!bootstrap){

    bsl = SuperLearner(Y = Y,
            X = Xmatb, 
            newX = newXmatb,
            family = family,
            SL.library = SL.library,
            method = method,
            id = id,
            verbose = verbose,
            control = control,
            cvControl = cvControl,
            obsWeights = weights)

      res = list(
      weights = bsl$coef,
      predictions = bsl$SL.predict,
      effectX = newXmatb
    )
      
    # nrow(res$effectX) == length(res$predictions)
    return(res)

  } else {
    
    stop("bootstrap not yet implemented.")
    # bootstrap SuperLearner
    # get effects for every combination of covariates
    # install.packages('doMC')
    # require('doMC')  # for easy parallelization on a nice multicore machine
    # registerDoMC(cores = cores)
    # bootres = foreach(i=1:R, .combine=rbind, .errorhandling='remove') %dopar% {
    res = list()
    for( i in 1:R){
      d = sample(x = 1:length(Xmatb[,1]),
                 size = length(Xmatb[,1]),
                 replace=TRUE)

      bsl = SuperLearner(Y = Y[d],
              X = Xmatb[d,],
              newX = newXb,
              family = family,
              SL.library = SL.library,
              method = method,
              id = id,
              verbose = verbose,
              control = control,
              cvControl = cvControl,
              obsWeights = weights)

       res[[i]] = list(weights = bsl$coef, predictions = bsl$SL.predict)

    }
  return(list(boostrap_samples = res, effectX = newXb))
  }
}

# TODO: process and format results, display for the user.




