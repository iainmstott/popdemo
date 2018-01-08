################################################################################
#' Determine reducibility of a matrix
#'
#' @description
#' Determine whether a matrix is irreducible or reducible
#'
#' @param A a square, non-negative numeric matrix of any dimension.
#'
#' @details 
#' \code{isIrreducible} works on the premise that a matrix \strong{A} 
#' is irreducible if and only if (\strong{I}+\strong{A})^(s-1) is positive, 
#' where \strong{I} is the identity matrix of the same dimension as \strong{A} 
#' and s is the dimension of \strong{A} (Caswell 2001).
#'
#' @return 
#' \code{TRUE} (for an irreducible matrix) or \code{FALSE} (for a reducible 
#' matrix).
#'
#' @references
#' Caswell (2001) matrix Population Models, 2nd. ed. Sinauer.
#'
#' @family PerronFrobeniusDiagnostics
#'
#' @examples
#'   # Create a 3x3 irreducible PPM
#'   ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#'
#'   # Diagnose reducibility
#'   isIrreducible(A)
#' 
#'   # Create a 3x3 reducible PPM
#'   B<-A; B[3,2] <- 0; B
#'
#'   # Diagnose reducibility
#'   isIrreducible(B)
#'
#' @concept
#' reducibility reducible irreducible Perron Frobenius
#'
#' @export isIrreducible
#' @importFrom expm "%^%"
#'
isIrreducible <-
function(A){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
I<-diag(order)
IplusA<-I+A
powermatrix<-IplusA%^%(order-1)
minval<-min(powermatrix)
if(minval>0){
    return(TRUE)
}
else{
    return(FALSE)
}}
