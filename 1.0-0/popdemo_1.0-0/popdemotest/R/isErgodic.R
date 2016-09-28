################################################################################
#' Determine ergodicity of a matrix
#'
#' @description
#' Determine whether a matrix is ergodic or nonergodic.
#'
#' @param A a square, non-negative numeric matrix of any dimension.
#'
#' @param digits the number of digits that the dominant left eigenvector should 
#' be rounded to.
#'
#' @param return.eigvec (optional) logical argument determining whether or not 
#' the dominant left eigenvector should be returned.
#'
#' @details 
#' \code{isErgodic} works on the premise that a matrix is ergodic if 
#' and only if the dominant left eigenvector (the reproductive value vector) of 
#' the matrix is positive (Stott \emph{et al}. 2010).\cr\cr
#' In rare cases, \R may calculate that the dominant left eigenvector of a 
#' nonergodic matrix contains very small entries that are approximate to (but 
#' not equal to) zero. Rounding the dominant eigenvector using \code{digits} 
#' prevents mistakes.\cr\cr
#'
#' @return 
#' If \code{return.eigvec=FALSE}, either \code{TRUE} (for an ergodic matrix) or 
#' \code{FALSE} (for a nonergodic matrix).\cr\cr
#' If \code{return.eigvec=TRUE}, a list containing elements:\cr
#' \item{\code{ergodic}}{\code{TRUE} or \code{FALSE}, as above.}
#' \item{\code{eigvec}}{the dominant left eigenvector of \code{A}.}
#'
#' @references
#' \href{http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2010.00032.x/full}{
#' Stott \emph{et al}. (2010) Methods. Ecol. Evol., 1, 242-252.}
#' \cr\cr
#'
#' @family Perron-Frobenius diagnostics
#'
#' @examples
#'   # Create a 3x3 ergodic PPM
#'   ( A <- matrix(c(0,0,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#'
#'   # Diagnose ergodicity
#'   isErgodic(A)
#'
#'   # Create a 3x3 nonergodic PPM
#'   B<-A; B[3,2] <- 0; B
#'
#'   # Diagnose ergodicity and return left eigenvector
#'   isErgodic(B, return.eigvec=TRUE)
#'
#' @concept
#' ergodicity ergodic nonergodic Perron Frobenius weak strong
#'
#' @export
#'
isErgodic <-
function(A,digits=12,return.eigvec=FALSE){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
leigs<-eigen(t(A))
lmax<-which.max(Re(leigs$values))
v<-leigs$vectors[,lmax]
Rev<-abs(Re(v))
Rev<-round(Rev,digits)
if(min(Rev)>0) ans<-TRUE else(ans<-FALSE)
if(min(Im(v))>0) leigvec<-v else(leigvec<-Rev)
if(return.eigvec){
    return(list(ergodic=ans,eigvec=leigvec))
}
else{
    return(ans)
}}
