################################################################################
#' Calculate elasticity matrix
#'
#' @description
#' Calculate the elasticity matrix for a specified population matrix projection 
#' model using eigenvectors.
#'
#' @param A a square, non-negative numeric matrix of any dimension, 
#' or a CompadreMat object (see RCompadre package).
#' @param eval the eigenvalue to evaluate. Default is \code{eval="max"}, which 
#' evaluates the dominant eigenvalue (the eigenvalue with largest REAL value: 
#' for imprimitive or reducible matrices this may not be the first eigenvalue). 
#' Otherwise, specifying e.g. \code{eval=2} will evaluate elasticity of the 
#' eigenvalue with second-largest modulus.

#' @details 
#' \code{elas} uses the eigenvectors of \code{A} to calculate the elasticity 
#' matrix of the specified eigenvalue, see section 9.1 in Caswell (2001). 
#' Same method as \code{elasticity} in \code{popbio} but can also evaluate
#' subdominant eigenvalues.
#'
#' @return 
#' A numeric (real or complex) matrix of equal dimension to \code{A}.
#'
#' @references
#' \itemize{
#'  \item Caswell (2001) Matrix Population Models 2nd ed. Sinauer.
#' }
#'
#' @family PerturbationAnalyses
#'
#' @examples
#'   # Create a 3x3 PPM
#'   ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#'
#'   # Calculate sensitivities of dominant eigenvalue
#'   elas(A)
#*
#'   # Calculate sensitivities of first subdominant eigenvalue,
#'   # only for observed transitions
#'   elas(A, eval=2)
#'
#' @concept 
#' perturbation sensitivity elasticity 
#'
#' @export elas
#' @importClassesFrom RCompadre CompadreMat
#' @importFrom RCompadre matA
#'
elas <-
function(A,eval="max"){
    if(class(A %in% "CompadreMat")){
        A <- matA(A)
    }
    if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
    order<-dim(A)[1]
    if(eval=="max"){
        val<-which.max(abs(Re(eigen(A)$values)))
    }
    else{
        val<-eval
    }
    reigs<-eigen(A)
    leigs<-eigen(t(A))
    lambda<-reigs$values[val]
    w<-as.matrix(reigs$vectors[,val])
    v<-as.matrix(leigs$vectors[,val])
    S<-(Conj(v)%*%t(w))/as.vector(Conj(t(v))%*%w)
    E<-(1/lambda)*S*A
    if(max(Im(E))>0) return(E) else(return(Re(E)))
}
