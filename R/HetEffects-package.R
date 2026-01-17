#' HetEffects: Ensemble Methods for Heterogeneous Treatment Effects
#'
#' The HetEffects package implements ensemble methods for estimating heterogeneous
#' treatment effects, extending Grimmer, Messing & Westwood (2017) with modern
#' causal machine learning methods.
#'
#' @name HetEffects-package
#' @aliases HetEffects
#' @keywords internal
#' @importFrom stats as.formula binomial complete.cases gaussian model.frame
#'   model.matrix model.response na.omit predict qnorm quantile sd setNames terms var
#' @importFrom graphics barplot hist lines
#' @importFrom grDevices rgb
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{het_ensemble}}: Fit an ensemble model for heterogeneous effects
#'   \item \code{\link{ate}}: Compute Average Treatment Effect
#'   \item \code{\link{cate}}: Compute Conditional Average Treatment Effect
#'   \item \code{\link{mcate}}: Compute Marginal Conditional Average Treatment Effect
#'   \item \code{\link{mate}}: Compute Marginal Average Treatment Effect
#' }
#'
#' @section Learners:
#' The package supports multiple machine learning methods:
#' \itemize{
#'   \item Core (glmnet): lasso, elastic_net_0.75, elastic_net_0.5, elastic_net_0.25
#'   \item Original paper: bayesglm, bart, randomforest, krls, svm
#'   \item Modern HTE: causal_forest (grf), xgboost
#' }
#'
#' @section Inference:
#' \itemize{
#'   \item \code{\link{analytic_se}}: Analytic standard errors via influence functions
#'   \item \code{\link{bootstrap_effects}}: Bootstrap inference
#'   \item \code{\link{confint.het_ensemble}}: Confidence intervals
#' }
#'
#' @section Observational Studies:
#' \itemize{
#'   \item \code{\link{estimate_propensity}}: Estimate propensity scores
#'   \item IPW weighting via \code{propensity} argument in \code{het_ensemble}
#' }
#'
#' @references
#' Grimmer, J., Messing, S., & Westwood, S. J. (2017). Estimating Heterogeneous
#' Treatment Effects and the Effects of Heterogeneous Treatments with Ensemble Methods.
#' \emph{Political Analysis}, 25(4), 413-434.
#'
"_PACKAGE"
