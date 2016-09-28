################################################################################
#' Calculate maximal amplification
#'
#' @description
#' Calculate maximal amplification for a population matrix projection model.
#'
#' @param A a square, primitive, non-negative numeric matrix of any dimension
#'
#' @param vector (optional) a numeric vector or one-column matrix describing 
#' the age/stage distribution ('demographic structure') used to calculate a 
#' 'case-specific' maximal amplification.
#'
#' @param return.N (optional) if \code{TRUE}, returns population size at the 
#' point of maximal amplification (including effects of asymptotic growth and 
#' initial population size), alongside standardised maximal amplification.
#'
#' @param return.t (optional) if \code{TRUE}, returns the time at which maximal 
#' amplification occurs in the population projection.
#'
#' @param return.stage (optional) if \code{TRUE} and \code{vector="n"}, returns 
#' the stage that achieves the bound on maximal amplification.
#'
#' @param conv.iterations the maximum number of iterations allowed when calulating 
#' convergence time (see details). Please see \code{iterations} in 
#' \code{\link{convt}}.
#'
#' @param conv.accuracy the accuracy of convergence (see details). Please see 
#' \code{accuracy} in \code{\link{convt}}.
#'
#' @details 
#' \code{maxamp} returns a standardised measure of maximal amplification, 
#' discounting the effects of both initial population size and asymoptotic growth 
#' (Stott et al. 2011).\cr\cr  
#' If \code{vector} is not specified then the bound on maximal amplification (the 
#' largest maximal amplification that may be achieved) is returned, otherwise a 
#' 'case-specific' maximal amplification for the specified matrix and demographic 
#' structure is calculated. Note that not all demographic structures will yield a 
#' maximal amplification: if the model does not amplify then an error is returned.\cr\cr
#' Setting \code{return.N=T}, \code{return.t=T} and \code{return.stage=T} results in 
#' the function returning realised population size at maximal amplification 
#' (including the effects of asymptotic growth and initial population size), the 
#' time at which maximal amplification occurs and (if \code{vector="n"}), 
#' the stage-bias that results in the bound on maximal amplification, respectively.
#' NOTE that \code{N} is not indicative of maximum possible population size for a 
#' non-standardised model: merely the population size at the point of maximal 
#' amplification (i.e. largest positive deviation from lambda-max).\cr\cr
#' \code{max.amp} uses a simulation technique, using \code{\link{project}} to project 
#' the dynamics of the model before evaluating maximum projected density over all t. 
#' \code{conv.accuracy} and \code{conv.iterations} are passed to 
#' \code{\link{convt}}, which is used to find the point of model convergence 
#' in order to ensure maximal amplification is correctly captured in model projection.\cr\cr
#' \code{maxamp} will not work for imprimitive or reducible matrices.
#'
#' @return 
#' If \code{vector="n"}, the bound on maximal amplification of \code{A}.\cr
#' If \code{vector} is specified, the case-specific maximal amplification of the model.\cr
#' If \code{return.N=TRUE}, \code{return.t=TRUE} and/or \code{return.stage=TRUE},
#' a list with possible components:\cr
#' \describe{
#' \item{maxamp}{the bound on or case-specific maximal amplification}
#' \item{N}{the population size at the point of maximal amplification, including the 
#' effects of initial population size and asymptotic growth. NOTE that \code{N} is not 
#' indicative of maximum possible population size for a non-standardised model:
#' merely the population size at the point of maximal amplification (i.e. largest 
#' positive deviation from lambda-max).}
#' \item{t}{the projection interval at which maximal amplification is achieved.}
#' \item{stage}{(only if \code{vector="n"}), the stage that achieves the bound on 
#' maximal amplification.}
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
#'   # Create an initial stage structure
#'   ( initial <- c(1,3,2) )
#'
#'   # Calculate the bound on maximal amplification of A
#'   maxamp(A)
#'
#'   # Calculate the bound on maximal amplification of A and 
#'   # return the stage that achieves it
#'   maxamp(A, return.stage=TRUE)
#'
#'   # Calculate case-specific maximal amplification of A
#'   # and initial
#'   maxamp(A, vector=initial)
#'
#'   # Calculate case-specific maximal amplification of A
#'   # and initial and return realised population size and the 
#'   # time at which it is achieved
#'   maxamp(A, vector=initial, return.N=TRUE, return.t=TRUE)
#'
#' @concept 
#' transient dynamics amplification unstable instability
#'
#' @export
#'
maxamp <-
function(A,vector="n",return.N=FALSE,return.t=FALSE,return.stage=FALSE,conv.iterations=1e+5,conv.accuracy=1e-5){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
M<-A
eigvals<-eigen(M)$values
lmax<-which.max(Re(eigvals))
lambda<-Re(eigvals[lmax])
A<-M/lambda
if(vector[1]=="n"){
    maxtime<-max(convt(A,accuracy=conv.accuracy,iterations=conv.iterations))
    projection<-project(A,time=maxtime)
    maxN<-numeric(order)
    times<-numeric(order)
    for(i in 1:order){
        maxN[i]<-max(projection[,i])
        times[i]<-which.max(projection[,i])-1
    }
    rhomax<-max(maxN)
    t<-times[which.max(maxN)]
    stage<-which.max(maxN)
    if(return.N){
        Nt<-rhomax*lambda^t
        if(all(return.t,return.stage)) return(list(maxamp=rhomax,N=Nt,t=t,stage=stage))
        if(all(return.t,!return.stage)) return(list(maxamp=rhomax,N=Nt,t=t))
        if(all(!return.t,return.stage)) return(list(maxamp=rhomax,N=Nt,stage=stage))
        if(!all(return.t,return.stage)) return(list(maxamp=rhomax,N=Nt))
    }
    else{
        if(all(return.t,return.stage)) return(list(maxamp=rhomax,t=t,stage=stage))
        if(all(return.t,!return.stage)) return(list(maxamp=rhomax,t=t))
        if(all(!return.t,return.stage)) return(list(maxamp=rhomax,stage=stage))
        if(!all(return.t,return.stage)) return(rhomax)
    }
}
else{
    n0<-vector
    vector<-n0/sum(n0)
    maxtime<-convt(A,vector=vector,accuracy=conv.accuracy,iterations=conv.iterations)
    projection<-project(A,vector=vector,time=maxtime)
    rhomax<-max(projection)
    t<-which.max(projection)-1
    if(rhomax>1){
        if(return.N){
            Nt<-rhomax*sum(n0)*lambda^t
            if(return.t) return(list(maxamp=rhomax,N=Nt,t=t)) else(return(list(maxamp=rhomax,N=Nt)))
        }
        else{
            if(return.t) return(list(maxamp=rhomax,t=t)) else(return(rhomax))
        }
    }
    else{
        stop("Model does not amplify. Cannot compute maximum amplification")
    }
}}
