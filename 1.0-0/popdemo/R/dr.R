################################################################################
#' Calculate damping ratio
#'
#' @description
#' Calculate the damping ratio of a given population matrix projection model.
#'
#' @param A a square, irreducible, non-negative numeric matrix of any dimension.
#'
#' @param return.time (optional) a logical argument determining whether an 
#' estimated convergence time should be returned.
#'
#' @param x (optional) the logarithm used in determining estimated time to 
#' convergence (see details).
#'
#' @details 
#' The damping ratio is calculated as the ratio of the dominant eigenvalue to 
#' the modulus of the largest subdominant eigenvalue. Time to convergence can 
#' be estmimated by calculating \code{log(dr)/log(x)}, which is the time taken 
#' for the dominant eigenvalue to become \code{x} times larger than the largest 
#' subdominant eigenvalue.
#'
#' @return 
#' If \code{return.time=FALSE}, the damping ratio of \code{A}.\cr
#' If \code{return.time=TRUE}, a list containing components:
#' \describe{
#' \item{dr}{the damping ratio of \code{A}}
#' \item{t}{the estimated time to convergence.}
#' }
#'
#' @references
#' Caswell (2001) Matrix Population Models 2nd. ed. Sinauer.\cr
#' Stott et al. (2010) Ecol. Lett., 14, 959-970.
#'
#' @family ConvergenceMeasures
#'
#' @examples
#'   # Create a 3x3 PPM
#'   A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
#'
#'   # Calculate damping ratio
#'   dr(A)
#'
#'   # Calculate damping ratio and time to convergence using a 
#'   # multiple of 10
#'   dr(A, return.time=TRUE, x=10)
#'
#' @concept 
#' converge convergence resilience stability
#'
#' @export
#'
dr <-
function(A,return.time=FALSE,x=10){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
if(!isIrreducible(A)) stop("Matrix is reducible")
if(!isPrimitive(A)) warning("Matrix is imprimitive")
eigvals<-eigen(A)$values
lmax<-which.max(Re(eigvals))
lambda<-Re(eigvals[lmax])
lambda2<-sort(Mod(eigvals),decreasing=TRUE)[2]
dr<-lambda/lambda2
if(return.time){
    t<-log(x)/log(dr)
    return(list(dr=dr,t=t))
}
else{
    return(dr)
}}
