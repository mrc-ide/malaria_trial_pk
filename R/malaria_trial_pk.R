#------------------------------------------------
#' @title Estimates pharmacokinetic parameters from trial data via MCMC
#'
#' @description Note, this package was forked from the MCMC package drjacoby, so
#'   there may be some names throughout that come from there.
#'
#' @docType package
#' @name malariatrialpk
NULL

#------------------------------------------------
#' @useDynLib malariatrialpk, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %>%
NULL

#------------------------------------------------
# unload DLL when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("malariatrialpk", libpath)  # nocov
}
