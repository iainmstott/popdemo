################################################################################
#' Deprecated functions in the popdemo package
#' 
#' @aliases
#' Cohen.cumulative Cohen.cumulative-deprecated 
#' convergence.time convergence.time-deprecated 
#' inertia.tfa inertia.tfa-deprecated 
#' inertia.tfamatrix inertia.tfamatrix-deprecated 
#' inertia.tfsens inertia.tfsens-deprecated 
#' inertia.tfsensmatrix inertia.tfsensmatrix-deprecated 
#' is.matrix_ergodic is.matrix_ergodic-deprecated 
#' is.matrix_irreducible is.matrix_irreducible-deprecated 
#' is.matrix_primitive is.matrix_primitive-deprecated 
#' Keyfitz.delta Keyfitz-delta-deprecated 
#' projection.distance projection.distance-deprecated 
#' tfa tfa-deprecated 
#' tfamatrix tfamatrix-deprecated 
#' tfsens tfsens-deprecated 
#' tfsensmatrix tfsensmatrix-deprecated 
#' minCS minCS-deprecated 
#' tf tf-deprecated 
#' reactivity reactivity-deprecated 
#' firststepatt firststepatt-deprecated 
#'
#' @rdname popdemo-deprecated
#' @name popdemo-deprecated
#' @docType package
#'
#' @description
#' Deprecated functions in the popdemo package
#'
#' @param ... Parameters to be passed to the new function versions
#'
#' @details
#' Many functions have become deprecated as of popdemo_1.0-0.
#' Most of these are straightforward renamings of functions. 
#' The old function names are provided for compatibility with older versions 
#' of 'popdemo', but may eventually be completely removed! 
#' Please update your code to use the new function names.\cr\cr
#' Most deprecated functions needed to be renamed because they included a period 
#' in the function name: the new function names don't use periods, which is a 
#' better approach for playing nicely with the S3 OO system (see Hadley Wickham's 
#' \href{http://adv-r.had.co.nz/OO-essentials.html}{OO field guide} for more
#' info). These are:\cr
#' \tabular{rl}{
#'   \code{Cohen.cumulative} \tab now called \code{CohenD}\cr
#'   \code{convergence.time} \tab now called \code{convt}\cr
#'   \code{inertia.tfa} \tab now called \code{tfa_inertia}\cr
#'   \code{inertia.tfamatrix} \tab now called \code{tfam_inertia}\cr
#'   \code{inertia.tfsens} \tab now called \code{tfs_inertia}\cr
#'   \code{inertia.tfsensmatrix} \tab now called \code{tfsm_inertia}\cr
#'   \code{is.matrix_ergodic} \tab now called \code{isErgodic}\cr
#'   \code{is.matrix_irreducible} \tab now called \code{isIrreducible}\cr
#'   \code{is.matrix_primitive} \tab now called \code{isPrimitive}\cr
#'   \code{Keyfitz.delta} \tab now called \code{KeyfitzD}\cr
#'   \code{projection.distance} \tab now called \code{projectionD}\cr
#' }
#' Some other functions have been renamed to keep consistency with new functions, 
#' and also to further avoid problems with S3 methods by making sure classes and 
#' functions don't have the same names:\cr
#' \tabular{rl}{
#'   \code{tfa} \tab now called \code{tfa_lambda}\cr
#'   \code{tfamatrix} \tab now called \code{tfam_lambda}\cr
#'   \code{tfsens} \tab now called \code{tfs_lambda}\cr
#'   \code{tfsensmatrix} \tab now called \code{tfsm_lambda}\cr
#' }
#' Some functions have been made internal:\cr
#' \tabular{rl}{
#'   \code{minCS} \tab now called \code{.minCS}\cr
#'   \code{tf} \tab now called \code{.tf}\cr
#' }
#' Two functions are deprecated because they have been merged into one:\cr
#' \tabular{rl}{
#'   \code{reactivity,firststepatt} \tab now handled by \code{\link{reac}}.\cr
#' }
#' Before, \code{reactivity} handled first-timestep amplification and 
#' \code{firststepatt} handled first-timestep attenuation. This is silly, because
#' a projection EITHER amplifies OR attenuates in the first timestep. Desptite 
#' the semantics, \code{reac} now deals with both amplification and attenuation
#' in the first timestep, everything that was calculable in the previous two
#' functions is also calculable in the one new function.
#' 
#' @export 
#'
Cohen.cumulative <- function(...) {
.Deprecated(msg=paste("'Cohen.cumulative' has been renamed and is deprecated.",
                      "The old name will continue to work for now, but may be removed in later versions.",
                      "Use the new name 'CohenD' instead. See help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
  CohenD(...)
}

convergence.time <- function(...) {
.Deprecated(msg=paste("'convergence.time' has been renamed and is deprecated.",
                      "The old name will continue to work for now, but may be removed in later versions.",
                      "Use the new name 'convt' instead. See help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
  convt(...)
}

inertia.tfa <- function(...) {
.Deprecated(msg=paste("'inertia.tfa' has been renamed and is deprecated.",
                      "The old name will continue to work for now, but may be removed in later versions.",
                      "Use the new name 'tfa_inertia' instead. See help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
  tfa_inertia(...)
}

inertia.tfamatrix <- function(...) {
.Deprecated(msg=paste("'inertia.tfamatrix' has been renamed and is deprecated.",
                      "The old name will continue to work for now, but may be removed in later versions.",
                      "Use the new name 'tfam_inertia' instead. See help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
  tfam_inertia(...)
}
inertia.tfsens <- function(...) {
.Deprecated(msg=paste("'inertia.tfsens' has been renamed and is deprecated.",
                      "The old name will continue to work for now, but may be removed in later versions.",
                      "Use the new name 'tfs_inertia' instead. See help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
  tfs_inertia(...)
}

inertia.tfsensmatrix <- function(...) {
.Deprecated(msg=paste("'inertia.tfsensmatrix' has been renamed and is deprecated.",
                      "The old name will continue to work for now, but may be removed in later versions.",
                      "Use the new name 'tfsm_inertia' instead. See help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
  tfsm_inertia(...)
}

is.matrix_ergodic <- function(...) {
.Deprecated(msg=paste("'is.matrix_ergodic' has been renamed and is deprecated.",
                      "The old name will continue to work for now, but may be removed in later versions.",
                      "Use the new name 'isErgodic' instead. See help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
  isErgodic(...)
}

is.matrix_irreducible <- function(...) {
.Deprecated(msg=paste("'is.matrix_irreducible' has been renamed and is deprecated.",
                      "The old name will continue to work for now, but may be removed in later versions.",
                      "Use the new name 'isIrreducible' instead. See help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
  isIrreducible(...)
}
is.matrix_primitive <- function(...) {
.Deprecated(msg=paste("'is.matrix_primitive' has been renamed and is deprecated.",
                      "The old name will continue to work for now, but may be removed in later versions.",
                      "Use the new name 'isPrimitive' instead. See help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
  isPrimitive(...)
}

Keyfitz.delta <- function(...) {
.Deprecated(msg=paste("'Keyfitz.delta' has been renamed and is deprecated.",
                      "The old name will continue to work for now, but may be removed in later versions.",
                      "Use the new name 'KeyfitzD instead. See help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
  KeyfitzD(...)
}

projection.distance <- function(...) {
.Deprecated(msg=paste("'projection.distance' has been renamed and is deprecated.",
                      "The old name will continue to work for now, but may be removed in later versions.",
                      "Use the new name 'projectionD' instead. See help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
  projectionD(...)
}

tfa <- function(...) {
.Deprecated(msg=paste("'tfa' has been renamed and is deprecated.",
                      "The old name will continue to work for now, but may be removed in later versions.",
                      "Use the new name 'tfa_lambda' instead. See help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
  tfa_lambda(...)
}

tfamatrix <- function(...) {
.Deprecated(msg=paste("'tfamatrix' has been renamed and is deprecated.",
                      "The old name will continue to work for now, but may be removed in later versions.",
                      "Use the new name 'tfam_lambda' instead. See help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
  tfam_lambda(...)
}

tfsens <- function(...) {
.Deprecated(msg=paste("'tfsens' has been renamed and is deprecated.",
                      "The old name will continue to work for now, but may be removed in later versions.",
                      "Use the new name 'tfs_lambda' instead. See help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
  tfs_lambda(...)
}

tfsensmatrix <- function(...) {
.Deprecated(msg=paste("'tfsensmatrix' has been renamed and is deprecated.",
                      "The old name will continue to work for now, but may be removed in later versions.",
                      "Use the new name 'tfsm_lambda' instead. See help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
  tfsm_lambda(...)
}

minCS <- function(...) {
.Deprecated(msg=paste("'minCS' has been made internal and is deprecated.",
                      "The old name will continue to work for now, but may be removed in later versions.",
                      "Use '.minCS' instead. See help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
  .minCS(...)
}

tf <- function(...) {
.Deprecated(msg=paste("'tf' has been made internal is deprecated.",
                      "The old name will continue to work for now, but may be removed in later versions.",
                      "Use the new name '.tf' instead. See help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
  .tf(...)
}

reactivity <-
function(...){
.Deprecated(msg=paste("'reactivity' is deprecated: it has been merged with 'firststepatt' into a single function.",
                      "Use the new function 'reac' instead. See help(\"reac\") and help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
reac(...)
}

firststepatt <-
function(...){
.Deprecated(msg=paste("'firststepatt' is deprecated: it has been merged with 'reactivity' into a single function.",
                      "Use the new function 'reac' instead. See help(\"reac\") and help(\"popdemo-deprecated\") for more info.",
                      sep="\n"))
reac(...)
}
