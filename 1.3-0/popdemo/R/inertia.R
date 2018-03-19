################################################################################
#' Calculate population inertia
#'
#' @description
#' Calculate population inertia for a population matrix projection model.
#'
#' @param A a square, primitive, irreducible, non-negative numeric matrix of any 
#' dimension
#' @param vector (optional) a numeric vector or one-column matrix describing 
#' the age/stage distribution ('demographic structure') used to calculate a 
#' 'case-specific' maximal amplification
#' @param bound (optional) specifies whether an upper or lower bound should be 
#' calculated (see details).
#' @param return.N (optional) if \code{TRUE}, returns population size for a 
#' specified \code{t} (including effects of asymptotic growth and initial 
#' population size), alongside standardised inertia.
#' @param t (optional) the projection interval at which \code{N} is to be 
#' calculated.  Calculation of \code{N} is only accurate for \code{t} where the 
#' model has converged (see details)
#'
#' @details 
#' A nonstable population, when it achieves asymptotic growth following transient 
#' dynamics, is a fixed ratio of the size of a population projected with the same
#' initial size but stable structure. \code{inertia} calculates the value of this
#' ratio (Koons et al. 2007)
#' 
#' If \code{vector="n"} then either \code{bound="upper"} or \code{bound="lower"}
#' must be specified, which calculate the upper or lower bound on population 
#' inertia (i.e. the largest and smallest values that inertia may take) 
#' respectively. Specifying \code{vector} overrides calculation of a bound, and 
#' will yield a 'case-specific' value for inertia.
#' 
#' \code{inertia} will not work with imprimitive or reducible matrices.
#'
#' @return 
#' If \code{vector="n"}, the upper bound on inertia of \code{A} if 
#' \code{bound="upper"} and the lower bound on inertia of \code{A} if 
#' \code{bound="lower"}.\cr
#' If \code{vector} is specified, the case-specific inertia of the model.\cr
#' If \code{return.N=TRUE} and \code{t} is specified, a list with components:
#' \describe{
#' \item{inertia}{the bound on or case-specific inertia}
#' \item{N}{the population size at specified \code{t}.}
#' }
#'
#' @references
#' \itemize{
#'  \item Koons et al. (2007) Ecology, 88, 2867-2867.
#'  \item Stott et al. (2011) Ecol. Lett., 14, 959-970.
#' }
#'
#' @family TransientIndices
#'
#' @seealso
#' Transfer function methods for inertia: \code{\link{inertia.tfa}}, 
#' \code{\link{inertia.tfamatrix}}, \code{\link{inertia.tfsens}}, 
#' \code{\link{inertia.tfsensmatrix}}
#'
#' @examples
#'   # Create a 3x3 PPM
#'   ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#'
#'   # Create an initial stage structure
#'   ( initial <- c(1,3,2) )
#'
#'   # Calculate the upper bound on inertia of A
#'   inertia(A,bound="upper")
#'
#'   # Calculate the lower bound on inertia of A
#'   inertia(A,bound="lower")
#'
#'   # Calculate case-specific inertia of A and initial
#'   inertia(A, vector=initial)
#'
#'   # Calculate case-specific inertia of A and initial and 
#'   # return realised population size at t=25
#'   inertia(A, vector=initial, return.N=TRUE, t=25)
#'
#' @concept 
#' transient amplification attenuation unstable instability stable equivalent ratio
#'
#' @export inertia
#'
inertia <-
function(A,vector="n",bound=NULL,return.N=FALSE,t=NULL){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
if(!isIrreducible(A)) stop("Matrix A is reducible")
if(!isPrimitive(A)) stop("Matrix A is imprimitive")
M<-A
reigs<-eigen(A)
leigs<-eigen(t(A))
lmax<-which.max(Re(reigs$values))
lambda<-Re(leigs$values[lmax])
A<-M/lambda
w<-as.matrix(abs(Re(reigs$vectors[,lmax])))
v<-as.matrix(abs(Re(leigs$vectors[,lmax])))
if(vector[1]=="n"){
    if(!any(bound=="upper",bound=="lower")) stop('Please specify bound="upper", bound="lower" or specify vector')
    if(bound=="upper"){
        rhoinfinity<-as.vector((max(v)*sum(w))/(t(v)%*%w))
        if(return.N){
            if(is.null(t)) stop("Please specify a value of t at which N is to be calculated")
            warning("Estimation of N will be  inaccurate for\n t where the model has not converged.")
            N<-rhoinfinity*lambda^t
            return(list(upper.inertia=rhoinfinity,N=N))
        }
        else{
            return(rhoinfinity)
        }
    }
    if(bound=="lower"){
        rhoinfinity<-as.vector((min(v)*sum(w))/(t(v)%*%w))
        if(return.N){
            if(is.null(t)) stop("Please specify a value of t at which N is to be calculated")
            warning("Estimation of N will be  inaccurate for\n t where the model has not converged.")
            N<-rhoinfinity*lambda^t
            return(list(lower.inertia=rhoinfinity,N=N))
        }
        else{
            return(rhoinfinity)
        }
    }
}
else{
    if(!is.null(bound)) warning("Specification of vector overrides calculation of bound")
    n0<-vector
    vector<-n0/sum(n0)
    Pinfinity<-as.vector((t(v)%*%vector*sum(w))/(t(v)%*%w))
    if(return.N){
        if(is.null(t)) stop("Please specify a value of t at which N is to be calculated")
        warning("Estimation of N will be  inaccurate for\n t where the model has not converged.")
        N<-Pinfinity*sum(n0)*lambda^t
        return(list(inertia=Pinfinity,N=N))
    }
    else{
        return(Pinfinity)
    }
}}

