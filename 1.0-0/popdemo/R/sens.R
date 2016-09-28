################################################################################
#' Calculate sensitivity matrix
#'
#' @description
#' Calculate the sensitivity matrix for a population matrix projection model 
#' using eigenvectors.
#'
#' @param A a square, non-negative numeric matrix of any dimension
#'
#' @param eval the eigenvalue to evaluate. Default is \code{eval="max"}, which 
#' evaluates the dominant eigenvalue (the eigenvalue with largest REAL value: 
#' for imprimitive or reducible matrices this may not be the first eigenvalue). 
#' Otherwise, specifying e.g. \code{eval=2} will evaluate sensitivity of the 
#' eigenvalue with second-largest modulus.
#'
#' @param all (optional) if \code{FALSE}, then only sensitivity values for 
#' observed transitions (nonzero entries in \code{A}) are returned.
#'
#' @details 
#' \code{sens} uses the eigenvectors of \code{A} to calculate the sensitivity 
#' matrix of the specified eigenvalue, see section 9.1 in Caswell (2001). 
#' Same method as \code{sensitivity} in \code{popbio} but can also evaluate
#' subdominant eigenvalues.
#'
#' @return 
#' A numeric (real or complex) matrix of equal dimension to \code{A}.
#'
#' @references
#' Caswell (2001) Matrix Population Models 2nd ed. Sinauer.
#'
#' @family PerturbationAnalyses
#'
#' @examples
#'   # Create a 3x3 PPM
#'   ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#'
#'   # Calculate sensitivities of dominant eigenvalue
#'   sens(A)
#*
#'   # Calculate sensitivities of first subdominant eigenvalue,
#'   # only for observed transitions
#'   sens(A, eval=2, all=FALSE)
#'
#' @concept 
#' perturbation sensitivity elasticity 
#'
#' @export
#'
sens <-
function(A,eval="max",all=FALSE){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
reigs<-eigen(A)
leigs<-eigen(t(A))
if(eval=="max"){
    val<-which.max(Re(reigs$values))
}
else{
    val<-eval
}
w<-as.matrix(reigs$vectors[,val])
v<-as.matrix(leigs$vectors[,val])
S<-(Conj(v)%*%t(w))/as.vector(Conj(t(v))%*%w)
if(!all) S[A==0]<-0
if(max(Im(S))>0) return(S) else(return(Re(S)))
}

