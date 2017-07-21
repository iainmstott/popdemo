################################################################################
#' Determine primitivity of a matrix
#'
#' @description
#' Determine whether a matrix is primitive or imprimitive
#'
#' @param A a square, non-negative numeric matrix of any dimension.
#'
#' @details 
#' \code{isPrimitive} works on the premise that a matrix \strong{A} is 
#' primitive if \strong{A}^(s^2-(2*s)+2) is positive, where s is the dimension 
#' of \strong{A} (Caswell 2001).
#'
#' @return 
#' \code{TRUE} (for an primitive matrix) or \code{FALSE} (for an imprimitive 
#' matrix).
#'
#' @references
#' Caswell (2001) matrix Population Models, 2nd. ed. Sinauer.
#'
#' @family PerronFrobeniusDiagnostics
#'
#' @examples
#'   # Create a 3x3 primitive PPM
#'   ( A <- matrix(c(0,1,2,0.5,0,0,0,0.6,0), byrow=TRUE, ncol=3) )
#'
#'   # Diagnose primitivity
#'   isPrimitive(A)
#'
#'   # Create a 3x3 imprimitive PPM
#'   B<-A; B[1,2] <- 0; B
#'
#'   # Diagnose primitivity
#'   isPrimitive(B)
#'
#' @concept
#' primitivity primitive imprimitive Perron Frobenius
#'
#' @export
#' @importFrom expm "%^%"
#'
isPrimitive <-
function(A){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
powermatrix<-A%^%((order^2)-(2*order)+2)
minval<-min(powermatrix)
if(minval>0){
    return(TRUE)
}
else{
    return(FALSE)
}}
