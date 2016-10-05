################################################################################
#' Calculate maximal attenuation
#'
#' @description
#' Calculate maximal attenuation for a population matrix projection model.
#'
#' @param A a square, primitive, non-negative numeric matrix of any dimension
#'
#' @param vector (optional) a numeric vector or one-column matrix describing 
#' the age/stage distribution ('demographic structure') used to calculate a 
#' 'case-specific' maximal attenuation
#'
#' @param return.N (optional) if \code{TRUE}, returns population size at the 
#' point of maximal attenuation (including effects of asymptotic growth and 
#' initial population size), alongside standardised maximal attenuation.
#'
#' @param return.t (optional) if \code{TRUE}, returns the time at which maximal 
#' attenuation occurs in the population projection.
#'
#' @param return.stage (optional) if \code{TRUE} and \code{vector="n"}, returns 
#' the stage that achieves the bound on maximal attenuation.
#'
#' @param conv.iterations the maximum number of iterations allowed when calulating 
#' convergence time (see details). Please see \code{iterations} in 
#' \code{\link{convt}}.
#'
#' @param conv.accuracy the accuracy of convergence (see details). Please see 
#' \code{accuracy} in \code{\link{convt}}.
#'
#' @details 
#' \code{maxatt} returns a standardised measure of maximal attenuation, 
#' discounting the effects of both initial population size and asymoptotic growth 
#' (Stott et al. 2011).\cr\cr  
#' If \code{vector} is not specified then the bound on maximal attenuation (the 
#' greatest maximal attenuation that may be achieved) is returned, otherwise a 
#' 'case-specific' maximal attenuation for the specified matrix and demographic 
#' structure is calculated. Note that not all demographic structures will yield a 
#' maximal attenuation: if the model does not amplify then an error is returned.\cr\cr
#' Setting \code{return.N=T}, \code{return.t=T} and \code{return.stage=T} results in 
#' the function returning realised population size at maximal attenuation 
#' (including the effects of asymptotic growth and initial population size), the 
#' time at which maximal attenuation occurs and (if \code{vector="n"}), 
#' the stage-bias that results in the bound on maximal attenuation, respectively.
#' NOTE that \code{N} is not indicative of minuium possible population size for a 
#' non-standardised model: merely the population size at the point of maximal 
#' attenuation (i.e. largest negative deviation from lambda-max).\cr\cr
#' \code{max.att} uses a simulation technique, using \code{\link{project}} to project 
#' the dynamics of the model before evaluating minimum projected density over all t. 
#' \code{conv.accuracy} and \code{conv.iterations} are passed to 
#' \code{\link{convt}}, which is used to find the point of model convergence 
#' in order to ensure maximal attenuation is correctly captured in model projection.\cr\cr
#' \code{maxatt} will not work for imprimitive or reducible matrices.
#'
#' @return 
#' If \code{vector="n"}, the bound on maximal attenuation of \code{A}.\cr
#' If \code{vector} is specified, the case-specific maximal attenuation of the model.\cr
#' If \code{return.N=TRUE}, \code{return.t=TRUE} and/or \code{return.stage=TRUE},
#' a list with possible components:\cr
#' \describe{
#' \item{maxatt}{the bound on or case-specific maximal attenuation}
#' \item{N}{the population size at the point of maximal attenuation, including the 
#' effects of initial population size and asymptotic growth. NOTE that \code{N} is not 
#' indicative of minimum possible population size for a non-standardised model:
#' merely the population size at the point of maximal attenuation (i.e. largest 
#' negative deviation from lambda-max).}
#' \item{t}{the projection interval at which maximal attenuation is achieved.}
#' \item{stage}{(only if \code{vector="n"}), the stage that achieves the bound on 
#' maximal attenuation.}
#' }
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
#'   ( initial <- c(3,1,1) )
#'
#'   # Calculate the bound on maximal attenuation of A
#'   maxatt(A)
#'
#'   # Calculate the bound on maximal attenuation of A and 
#'   # return the stage that achieves it
#'   maxatt(A, return.stage=TRUE)
#'
#'   # Calculate case-specific maximal attenuation of A
#'   # and initial
#'   maxatt(A, vector=initial)
#'
#'   # Calculate case-specific maximal attenuation of A
#'   # and initial and return realised population size and the 
#'   # time at which it is achieved
#'   maxatt(A, vector=initial, return.N=TRUE, return.t=TRUE)
#'
#' @concept 
#' transient dynamics attenuation unstable instability
#'
#' @export
#'
maxatt <-
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
    minN<-numeric(order)
    times<-numeric(order)
    for(i in 1:order){
        minN[i]<-min(projection[,i])
        times[i]<-which.min(projection[,i])-1
    }
    rhomin<-min(minN)
    t<-times[which.min(minN)]
    stage<-which.min(minN)
    if(return.N){
        Nt<-rhomin*lambda^t
        if(all(return.t,return.stage)) return(list(maxatt=rhomin,N=Nt,t=t,stage=stage))
        if(all(return.t,!return.stage)) return(list(maxatt=rhomin,N=Nt,t=t))
        if(all(!return.t,return.stage)) return(list(maxatt=rhomin,N=Nt,stage=stage))
        if(!all(return.t,return.stage)) return(list(maxatt=rhomin,N=Nt))
    }
    else{
        if(all(return.t,return.stage)) return(list(maxatt=rhomin,t=t,stage=stage))
        if(all(return.t,!return.stage)) return(list(maxatt=rhomin,t=t))
        if(all(!return.t,return.stage)) return(list(maxatt=rhomin,stage=stage))
        if(!all(return.t,return.stage)) return(rhomin)
    }
}
else{
    n0<-vector
    vector<-n0/sum(n0)
    maxtime<-convt(A,vector=vector,accuracy=conv.accuracy,iterations=conv.iterations)
    projection<-project(A,vector=vector,time=maxtime)
    rhomin<-min(projection)
    t<-which.min(projection)-1
    if(rhomin<1){
        if(return.N){
            Nt<-rhomin*sum(n0)*lambda^t
            if(return.t) return(list(maxatt=rhomin,N=Nt,t=t)) else(return(list(maxatt=rhomin,N=Nt)))
        }
        else{
            if(return.t) return(list(maxatt=rhomin,t=t)) else(return(rhomin))
        }
    }
    else{
        stop("Model does not attenuate.  Cannot compute maximum attenuation")
    }
}}

