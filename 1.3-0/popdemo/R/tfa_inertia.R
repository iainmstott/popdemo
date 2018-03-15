################################################################################
#' Transfer function Analysis
#'
#' @description
#' Transfer function analysis of inertia of a population matrix 
#' projection model using a specified perturbation structure.
#'
#' @param A a square, primitive, nonnegative numeric matrix of any dimension
#'
#' @param d,e numeric vectors that determine the perturbation structure 
#' (see details).
#'
#' @param vector (optional) a numeric vector or one-column matrix describing 
#' the age/stage distribution ('demographic structure') used to calculate the
#' transfer function of a 'case-specific' inertia
#'
#' @param bound (optional) specifies whether the transfer funciton of an upper 
#' or lower bound on inertia should be calculated (see details).
#'
#' @param prange a numeric vector giving the range of perturbation magnitude
#' (see details)
#'
#' @param digits specifies which values of lambda should be excluded from 
#' analysis to avoid a computationally singular system (see details).
#'
#' @details 
#' \code{tfa_inertia} calculates the transfer function of inertia of a 
#' population matrix projection model given a perturbation structure 
#' (specified using \code{d} and \code{e}), and a range of desired perturbation
#' magnitude (\code{prange}). Currently \code{tfa_inertia} can only work with 
#' rank-one, single-parameter perturbations (see Hodgson & Townley 2006).\cr\cr
#' If \code{vector="n"} then either \code{bound="upper"} or \code{bound="lower"}
#' must be specified, which calculate the transfer function of the upper or 
#' lower bound on population inertia (i.e. the largest and smallest values that 
#' inertia may take) respectively.  Specifying \code{vector} overrides 
#' calculation of a bound, and will yield a transfer function of a 'case-specific'
#' inertia.\cr\cr
#' The perturbation structure is determined by \code{d\%*\%t(e)}. Therefore, 
#' the rows to be perturbed are determined by \code{d} and the columns to be 
#' perturbed are determined by \code{e}. The specific values in d and e will 
#' determine the relative perturbation magnitude. So for example, if only entry
#' [3,2] of a 3 by 3 matrix is to be perturbed, then \code{d = c(0,0,1)} and 
#' \code{e = c(0,1,0)}. If entries [3,2] and [3,3] are to be perturbed with the 
#' magnitude of perturbation to [3,2] half that of [3,3] then \code{d = c(0,0,1)} 
#' and \code{e = c(0,0.5,1)}. \code{d} and \code{e} may also be expressed as 
#' numeric one-column matrices, e.g. \code{d = matrix(c(0,0,1), ncol=1)}, 
#' \code{e = matrix(c(0,0.5,1), ncol=1)}. See Hodgson et al. (2006) for more 
#' information on perturbation structures.\cr\cr
#' The perturbation magnitude is determined by \code{prange}, a numeric vector 
#' that gives the continuous range of perterbation magnitude to evaluate over. 
#' This is usually a sequence (e.g. \code{prange=seq(-0.1,0.1,0.001)}), but 
#' single transfer functions can be calculated using a single perturbation 
#' magnitude (e.g. \code{prange=-0.1}). Because of the nature of the equation
#' used to calculate the transfer function, \code{prange} is used to find a 
#' range of lambda from which the perturbation magnitudes are back-calculated, 
#' and matched to their orresponding inertia, so the output perturbation
#' magnitude \code{p} will match \code{prange} in length and range but not in
#' numerical value (see Stott et al. 2012 for more information).\cr\cr
#' \code{tfa_inertia} uses the resolvent matrix in its calculation, which
#' cannot be computed if any lambda in the equation are equal to the dominant
#' eigenvalue of \code{A}. \code{digits} specifies the values of lambda that 
#' should be excluded in order to avoid a computationally singular system. 
#' Any values equal to the dominant eigenvalue of \code{A} rounded to an 
#' accuracy of \code{digits} are excluded. \code{digits} should only need to 
#' be changed when the system is found to be computationally singular, in 
#' which case increasing \code{digits} should help to solve the problem.\cr\cr
#' \code{tfa_inertia} will not work for reducible matrices.\cr\cr
#' There is an S3 plotting method available (see \code{\link{plot.tfa}} 
#' and examples below)
#' 
#' @return 
#' A list containing numerical vectors:
#' \describe{
#' \item{p}{perturbation magnitudes}
#' \item{lambda}{dominant eigenvalues of perturbed matrices}
#' item{inertia}{inertias of perturbed matrices}
#' }
#' (Note that \code{p} will not be exactly the same as \code{prange} when 
#' \code{prange} is specified, as the code calculates p for a given lambda 
#' rather than the other way around, with \code{prange} only used to determine 
#' max, min and number of lambda values to evaluate.)
#'
#' @references
#' Stott et al. (2012) Methods Ecol. Evol., 3, 673-684.
#' Hodgson et al. (2006) J. Theor. Biol., 70, 214-224.
#'
#' @family TransferFunctionAnalyses
#' @family PerturbationAnalyses
#'
#' @seealso
#' S3 plotting method: \code{\link{plot.tfa}}\cr
#'
#' @examples
#'   # Create a 3x3 matrix
#'   ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#' 
#'   # Create an initial stage structure
#'   ( initial <- c(1,3,2) )
#' 
#'   # Calculate the transfer function of upper bound on inertia 
#'   # given a perturbation to A[3,2]
#'   ( transfer<-tfa_inertia(A, d=c(0,0,1), e=c(0,1,0), bound="upper",
#'                           prange=seq(-0.6,0.4,0.01)) )
#' 
#'   # Plot the transfer function using the S3 method (defaults to p 
#'   # and inertia in this case)
#'   plot(transfer)
#'
#'   # Plot inertia against lambda using the S3 method
#'   plot(transfer, xvar="lambda", yvar="inertia")
#' 
#'   # Calculate the transfer function of case-specific inertia 
#'   # given perturbation to A[3,2] and A[3,3] with perturbation 
#'   # to A[3,2] half that of A[3,3]
#'   ( transfer2<-tfa_inertia(A, d=c(0,0,1), e=c(0,0.5,1), vector=initial,
#'                            prange=seq(-0.6,0.4,0.01)) )
#' 
#'   # Plot inertia against p using the S3 method
#'   plot(transfer2)
#' 
#'   # Plot inertia against lambda without using the S3 method
#'   plot(transfer$inertia~transfer$lambda,type="l", 
#'        xlab=expression(lambda),ylab="inertia")
#' 
#' @concept 
#' KEYWORDS
#'
#' @export tfa_inertia
#'
tfa_inertia <-
function(A,d,e,vector="n",bound=NULL,prange,digits=1e-10){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!isIrreducible(A)) stop("Matrix A is reducible")
if(!isPrimitive(A)) stop("Matrix A is imprimitive")
eigvals<-eigen(A)$values
lmax<-which.max(Re(eigvals))
lambda<-Re(eigvals[lmax])
if(any(is.null(d),is.null(e))) stop("please specify a perturbation structure using d and e")
d<-as.matrix(d)
e<-as.matrix(e)
mineigvals<-eigen(A+(min(prange)*d%*%t(e)))$values
maxeigvals<-eigen(A+(max(prange)*d%*%t(e)))$values
minlmax<-which.max(Re(mineigvals))
minlambda<-Re(mineigvals[minlmax])
maxlmax<-which.max(Re(maxeigvals))
maxlambda<-Re(maxeigvals[maxlmax])
lambdarange<-seq(minlambda,maxlambda,(maxlambda-minlambda)/(length(prange)-1))
lambdarange<-lambdarange[round(lambdarange,-log10(digits))!=round(lambda,-log10(digits))]
pvals<-1/.tf(A,lambdarange,d,e)
I<-diag(order)
c<-as.matrix(rep(1,order))
bigd<-matrix(c(0,1),ncol=1)%x%d
bige<-matrix(c(1,0),ncol=1)%x%e
bigA2<-(matrix(c(1,0,0,1),ncol=2,byrow=TRUE)%x%A)+(matrix(c(0,1,0,0),ncol=2,byrow=TRUE)%x%I)
if(vector[1]=="n"){
    if(!any(bound=="upper",bound=="lower")) stop('Please specify bound="upper", bound="lower" or specify vector')
    vrange<-as.list(numeric(length(lambdarange)))
    mvrange<-numeric(length(lambdarange))
    vectorrange<-rep(list(as.matrix(numeric(order))),length(lambdarange))
    bigA1range<-as.list(numeric(length(lambdarange)))
    inertia<-numeric(length(lambdarange))
    if(bound=="upper"){
        for(i in 1:length(lambdarange)){
            vrange[[i]]<-abs(t(e)%*%solve((lambdarange[i]*I)-A))
            mvrange[i]<-which.max(vrange[[i]])
            vectorrange[[i]][mvrange[i]]<-1
            bigA1range[[i]]<-(matrix(c(1,0,0,1),ncol=2,byrow=TRUE)%x%A)+(matrix(c(0,1,0,0),ncol=2,byrow=TRUE)%x%(vectorrange[[i]]%*%t(c)))
            inertia[i]<-.tf(bigA1range[[i]],lambdarange[i],bigd,bige)/.tf(bigA2,lambdarange[i],bigd,bige)
        }
    }
    if(bound=="lower"){
        for(i in 1:length(lambdarange)){
            vrange[[i]]<-abs(t(e)%*%solve((lambdarange[i]*I)-A))
            mvrange[i]<-which.min(vrange[[i]])
            vectorrange[[i]][mvrange[i]]<-1
            bigA1range[[i]]<-(matrix(c(1,0,0,1),ncol=2,byrow=TRUE)%x%A)+(matrix(c(0,1,0,0),ncol=2,byrow=TRUE)%x%(vectorrange[[i]]%*%t(c)))
            inertia[i]<-.tf(bigA1range[[i]],lambdarange[i],bigd,bige)/.tf(bigA2,lambdarange[i],bigd,bige)
        }
    }
}
else{
    vector<-vector/sum(vector)
    bigA1<-(matrix(c(1,0,0,1),ncol=2,byrow=TRUE)%x%A)+(matrix(c(0,1,0,0),ncol=2,byrow=TRUE)%x%(vector%*%t(c)))
    inertia<-numeric(length(lambdarange))
    for(i in 1:length(lambdarange)){
        inertia[i]<-.tf(bigA1,lambdarange[i],bigd,bige)/.tf(bigA2,lambdarange[i],bigd,bige)
    }
}
final<-list(p=pvals,lambda=lambdarange,inertia=inertia)
class(final)<-c("tfa","list")
return(final)
}
