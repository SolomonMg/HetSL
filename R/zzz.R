#' @keywords internal
.onLoad <- function(libname, pkgname) {
  .init_learner_registry()
  invisible()
}
