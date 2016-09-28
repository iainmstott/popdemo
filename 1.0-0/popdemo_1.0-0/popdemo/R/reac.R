################################################################################
#' Calculate reactivity and first-timestep attenuation
#'
#' @description
#' Calculate reactivity (first-timestep amplification) and first-timestep 
#' attenuation for a population matrix projection model.
#'
#' @param A a square, non-negative numeric matrix of any dimension
#'
#' @param vector (optional) a numeric vector or one-column matrix describing 
#' the age/stage distribution used to calculate a 'case-specific' reactivity/
#' first-timestep attenuation
#'
#' @param return.N (optional) if \code{TRUE}, returns population size in the 
#' first time interval (including effects of asymptotic growth and initial 
#' population size), alongside standardised reactivity/first-timestep attenuation.
#'
#' @details 
#' \code{reac} returns a standardised measure of first-timestep 
#' amplification or attenuation, discounting the effects of both initial 
#' population size and asymoptotic growth (Stott et al. 2011).\cr\cr
#' If \code{vector}="n" then either \code{bound="upper"} or \code{bound="lower"}
#' must be specified, which calculate the upper or lower bound on first-timestep
#' amplification and attenuation (i.e. the largest and smallest values that 
#' reactivity and first-timestep attenuation may take) respectively.
#' Specifying \code{vector} overrides calculation of a bound, and will yield 
#' a 'case-specific' reactivity/first-timestep attenuation.\cr\cr
#' If \code{return.N=T} then the function also returns realised population size 
#' (including the effects of asymptotic growth and initial population size).\cr\cr
#' \code{reac} works with imprimitive and irreducible matrices, but 
#' returns a warning in these cases.\cr\cr
#' NOTE: \code{reac} replaces \code{reactivity} and \code{firststepatt} as of 
#' version 1.0-0. Although semantically 'reactivity' and 'first-timestep 
#' attenuation' are different (the former is an amplification in the first timestep
#' and the latter an attenuation in the first timestep), as a population matrix 
#' projection model EITHER amplifies OR attenuates in the first timestep, it made 
#' no sense to have two separate functions to calculate one thing 
#' (transient dynamics in the first timestep).
#'
#' @return 
#' If \code{vector="n"}, the upper bound on reactivity of \code{A} if 
#' \code{bound="upper"} and the lower bound on first-timestep attenuation of 
#' \code{A} if \code{bound="lower"}.\cr
#' If \code{vector} is specified, the 'case-specific' reactivity or first-timestep 
#' attenuation of the model.\cr\cr
#' If \code{return.N=TRUE}, a list with components:
#' \describe{
#' \item{reac}{the bound on or case-specific reactivity or first-timestep 
#' attenuation}
#' \item{N}{the population size at the first timestep, including the effects 
#' of initial population size and asymptotic growth.}
#'
#' @references
#' Neubert & Caswell (1997) Ecology, 78, 653-665.\cr
#' Stott et al. (2011) Ecol. Lett., 14, 959-970.\cr
#' Townley & Hodgson (2008) J. Appl. Ecol., 45, 1836-1839.
#'
#' @family TransientIndices
#'
#' @examples
#'   # Create a 3x3 PPM
#'   ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#'
#'   # Create initial stage structures
#'   ( initial1 <- c(1,3,2) )
#'   ( initial2 <- c(3,1,1) )
#'
#'   # Calculate the upper bound on reactivity of A
#'   reac(A, bound="upper")
#'
#'   # Calculate the lower bound on first-timestep attenuation of A
#'   reac(A, bound="lower")
#'
#'   # Calculate case-specific reactivity of A
#'   # when projected using specific demographic structure
#'   # that amplifies
#'   reac(A, vector=initial1)
#'
#'   # Calculate case-specific reactivity of A
#'   # and initial1 and return realised population size
#'   reac(A, vector=initial1, return.N=TRUE)
#'
#'   # Calculate case-specific first-timestep attenuation of A 
#'   # when projected using a specific demographic structure that¨
#'   #attenuates
#'   reac(A, vector=initial2)
#'
#'   # Calculate case-specific first-timestep attenuation of A 
#'   # and initial2 and return realised population size
#'   reac(A, vector=initial2, return.N=TRUE)#'
#'
#' @concept 
#' transient dynamics amplification unstable instability structure
#'
#' @export
#'
reac <-
function(A,vector="n",return.N=FALSE){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!is.matrix_irreducible(A)){
    warning("Matrix is reducible")
}
else{
    if(!is.matrix_primitive(A)) warning("Matrix is imprimitive")
}
M<-A
eigvals<-eigen(M)$values
lmax<-which.max(Re(eigvals))
lambda<-Re(eigvals[lmax])
A<-M/lambda
if(vector[1]=="n"){
    if(!any(bound=="upper",bound=="lower")) stop('Please specify bound="upper", bound="lower" or specify vector')
    if(bound=="upper"){
        reac<-norm(A)
        if(return.N){
            N<-reac*lambda
            return(list(reac=reac,N=N))
        }
        else{
            return(reac)
        }
    }
    if(bound=="lower"){
        reac<-.minCS(A)
        if(return.N){
            N<-reac*lambda
            return(list(reac=reac,N=N))
        }
        else{
            return(reac)
        }
    }
}
else{
    if(!is.null(bound)) warning("Specification of vector overrides calculation of bound")
    n0<-vector
    vector<-n0/sum(n0)
    reac<-sum(A%*%vector)
    if(return.N){
        N<-reac*sum(n0)*lambda
        return(list(reac=reac,N=N))
    }
    else{
        return(reac)
    }
}}
