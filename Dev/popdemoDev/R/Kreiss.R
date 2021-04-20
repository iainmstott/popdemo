################################################################################
#' Calculate Kreiss bounds
#'
#' @description
#' Calculate the upper or lower Kreiss bound for a population matrix projection 
#' model.
#'
#' @param A a square, irreducible, non-negative numeric matrix of any dimension, 
#' or a CompadreMat object (see RCompadre package).
#' @param bound (optional) specifies whether an upper or lower bound should be 
#' calculated.
#' @param return.r (optional) specifies whether the value of r at which the 
#' Kreiss bound is achieved should be returned (see details).
#' @param theta the value to which the Kriess bound is to be assessed relative 
#' to (see details).
#' @param rlimit the maximum value of r that may be reached before the code 
#' breaks (see details).
#' @param step1,step2 determine the iterative process in calculating the Kreiss 
#' bound (see details).
#'
#' @details 
#'  \code{Kreiss} by default returns a standardised Kreiss bound relative to both 
#' asymptotic growth/decline and initial population density (Townley & Hodgson 2008; 
#' Stott et al. 2011).  It uses an iterative process that evaluates a function of 
#' the resolvent of \code{A} over a range of values r where r>\code{theta}. This 
#' iterative process finds the maximum/minimum of the function for the upper/lower 
#' bounds respectively. The process is determined using \code{step1} and 
#' \code{step2}: in order to increase accuracy but keep computation time low, the 
#' function is evaluated forward in steps equal to \code{step1} until the 
#' maximum/minimum is passed and then backward in steps of \code{step2} to more 
#' accurately find the maximum/minimum itself. Therefore, \code{step1} should be 
#' larger than \code{step2}. The balance between both will determine computation 
#' time, whilst accuracy is determined almost solely by \code{step2}. The defaults 
#' should be sufficient for most matrices.
#' 
#' \code{theta} defaults to 1, which means the Kriess bound is assessed relative to 
#' both asymptotic growth and initial population size. Sometimes, the maximum/minimum 
#' of the function occurs at r-->\code{theta}, in which case r is equal to 
#' \code{theta+step2}. Setting \code{return.r=TRUE} tells the function to return the 
#' value of r where the maximum/minimum occurs alongside the value of the Kreiss bound. 
#' r may not exceed \code{rlimit}.
#' 
#' \code{Kreiss} will not work with reducible matrices, and returns a warning for 
#' imprimitive matrices.
#'
#' @return 
#' The upper or lower Kreiss bound of \code{A}.\cr
#' If \code{return.r=TRUE}, a list with components:
#' \describe{
#' \item{bound}{the upper or lower Kriess bound}
#' \item{r}{the value of r at which the function is minimised/maximised.}
#' }
#'
#' @references
#' \itemize{
#'  \item Stott et al. (2011) Ecol. Lett., 14, 959-970.
#'  \item Townley & Hodgson (2008) J. Appl. Ecol., 45, 1836-1839.
#' }
#'
#' @family TransientIndices
#'
#' @examples
#'   # Create a 3x3 PPM
#'   ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#'
#'   # Calculate the upper Kreiss bound of A
#'   Kreiss(A, bound="upper")
#'
#'   # Calculate the lower Kreiss bound of A
#'   Kreiss(A, bound="lower")
#'
#'   # Calculate the upper Kreiss bound of A and return 
#'   # the value of r at which the function is maximised
#'   Kreiss(A, bound="upper", return.r=TRUE)
#'
#' @concept 
#' transient amplification attenuation systems control unstable instability
#'
#' @export Kreiss
#' @importClassesFrom RCompadre CompadreMat
#' @importFrom RCompadre matA
#'
Kreiss <-
function(A,bound=NULL,return.r=FALSE,theta=1,rlimit=100,step1=1e-3,step2=1e-6){
    if(class(A %in% "CompadreMat")){
        A <- matA(A)
    }
    if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
    order<-dim(A)[1]
    if(!isIrreducible(A)) stop("Matrix A is reducible")
    if(!isPrimitive(A)) warning("Matrix A is imprimitive")
    M<-A
    eigvals<-eigen(M)$values
    lmax<-which.max(Re(eigvals))
    lambda<-Re(eigvals[lmax])
    A<-M/lambda
    I<-diag(order)
    if(is.null(bound)) stop('Please specify either bound="upper or bound="lower"')
    if(bound=="upper"){
        r1<-theta+step1
        K1<-(r1-theta)*(norm(solve((r1*I)-A)))
        Kreissbound1<-K1
        while(K1>=Kreissbound1 & r1<rlimit){
            Kreissbound1<-K1
            r1<-r1+step1
            K1<-(r1-theta)*(norm(solve((r1*I)-A)))
        }
        if(!r1<rlimit) stop("Maximum not reached and rlimit exceeded")
        r2<-r1
        K2<-(r2-theta)*(norm(solve((r2*I)-A)))
        Kreissbound2<-K1
        while(K2>=Kreissbound2 & r2>=(theta+(2*step2))){
            Kreissbound2<-K2
            r2<-r2-step2
            K2<-(r2-theta)*(norm(solve((r2*I)-A)))
        }
        rstar<-r2+step2
        Kstar<-Kreissbound2
        r.1<-r2
        K.1<-K2
        if(r2<(theta+(2*step2))){
            Kreissbound<-K.1
            r<-r.1
        }
        else{
            Kreissbound<-Kstar
            r<-rstar
        }
        if(return.r) return(list(Kreissbound=Kreissbound,r=r)) else(return(Kreissbound))
    }
    if(bound=="lower"){
        r1<-theta+step1
        K1<-(r1-theta)*(.minCS(solve((r1*I)-A)))
        Kreissbound1<-K1
        while(K1<=Kreissbound1 & r1<rlimit){
            Kreissbound1<-K1
            r1<-r1+step1
            K1<-(r1-theta)*(.minCS(solve((r1*I)-A)))
        }
        if(!r1<rlimit) stop("Minimum not reached and rlimit exceeded")
        r2<-r1
        K2<-(r2-theta)*(.minCS(solve((r2*I)-A)))
        Kreissbound2<-K1
        while(K2<=Kreissbound2 & r2>=(theta+(2*step2))){
            Kreissbound2<-K2
            r2<-r2-step2
            K2<-(r2-theta)*(.minCS(solve((r2*I)-A)))
        }
        rstar<-r2+step2
        Kstar<-Kreissbound2
        r.1<-r2
        K.1<-K2
        if(r2<(theta+(2*step2))){
            Kreissbound<-K.1
            r<-r.1
        }
        else{
            Kreissbound<-Kstar
            r<-rstar
        }
        if(return.r) return(list(Kreissbound=Kreissbound,r=r)) else (return(Kreissbound))
    }
}
