################################################################################
#' Transfer function Analysis
#'
#' @description
#' Transfer function analysis of inertia of a population matrix 
#' projection model for all matrix elements.
#'
#' @param A a square, primitive, nonnegative numeric matrix of any dimension
#'
#' @param vector (optional) a numeric vector or one-column matrix describing 
#' the age/stage distribution ('demographic structure') used to calculate the
#' transfer function of a 'case-specific' inertia
#'
#' @param bound (optional) specifies whether the transfer funciton of an upper 
#' or lower bound on inertia should be calculated (see details).
#'
#' @param elementtype (optional) a character matrix of the same dimension as 
#' \code{A} describing the structure of \code{A}: \code{"P"} denotes elements 
#' bounded between 0 and 1, i.e. survival, growth, regression; \code{"F"} 
#' denotes elements not bounded at 1, i.e. fecundity, fission; \code{NA} 
#' denotes absent elements (see details).
#'
#' @param Flim,Plim the perturbation ranges for \code{"F"} and \code{"P"} 
#' elements, expressed as a proportion of their magnitude  (see details).
#'
#' @param plength the desired length of the perturbation ranges.
#'
#' @param digits specifies which values of lambda should be excluded from 
#' analysis to avoid a computationally singular system (see details).
#'
#' @details 
#' \code{tfam_inertia} calculates an array of transfer functions 
#' of population inertia. A separate transfer function for each nonzero 
#' element of \code{A} is calculated (each element perturbed independently of 
#' the others). The function is most useful for use with the S3 method 
#' \code{\link{plot.tfam}} to visualise how perturbations affect the 
#' life cycle transitions, and easily compare the (nonlinear) effect of 
#' perturbation to different transitions on the dominant eigenvalue.\cr\cr
#' The sizes of the perturbations are determined by \code{elementtype}, 
#' \code{Flim}, \code{Plim} and \code{plength}. \code{elementtype} gives the 
#' type of each element, specifying whether perturbations should be 
#' bounded at 1 (\code{elementtype = "P"}) or not (\code{elementtype = "F"}). 
#' If \code{elementtype} is not directly specified, the function assigns its 
#' own types, with those in the first row attributed \code{"F"}, and elsewhere 
#' in the matrix attributed \code{"F"} if the value of the element >1 and 
#' \code{"P"} if the value of the element is <=1. \code{Flim} and \code{Plim} 
#' determine the desired perturbation magnitude, expressed as a proportion of 
#' the magnitude of the elements of \code{A}, whilst plength determines the 
#' length of the perturbation vector.  For example, if an "F" element is equal 
#' to 0.5, \code{Flim=c(-1,10)} and \code{plength=100} then the perturbation 
#' to that element is \code{seq(-1*0.5,10*0.5,100-1)}. The process is the same 
#' for \code{"P"} elements, except that these are truncated to a maximum value 
#' of 1 (growth/survival elements cannot be greater than 1). Both \code{"F"} 
#' and \code{"P"} elements are truncated to a minimum value of 0.\cr\cr
#' \code{tfam_inertia} uses \code{\link{tfa_inertia}} to calculate
#' transfer functions. \code{digits} is passed to \code{tfa_inertia} to
#' prevent the problem of singular matrices (see details in
#' \code{\link{tfa_inertia}}).\cr\cr
#' \code{tfam_inertia} will not work for reducible matrices.
#' 
#' @return 
#' A list containing numerical arrays: 
#' \describe{
#' \item{p}{perturbation magnitudes}
#' \item{lambda}{dominant eigenvalues of perturbed matrices}
#' \item{inertia}{inertias of perturbed matrices}
#' }
#' The first and second dimensions of the arrays are equivalent to the 
#' first and second dimensions of \code{A}. The third dimension of the 
#' arrays are the vectors returned by \code{tfa_inertia}. e.g.
#' $inertia[3,2,] selects the inertia values for the transfer function of
#' element [3,2] of the matrix.
#'
#' @references
#' Stott et al. (2012) Methods Ecol. Evol., 3, 673-684.
#' Hodgson et al. (2006) J. Theor. Biol., 70, 214-224.
#'
#' @family TransferFunctionAnalyses
#' @family PerturbationAnalyses
#'
#' @seealso
#' S3 plotting method:\cr
#' \code{\link{plot.tfam}}\cr
#'
#' @examples
#'   # Create a 3x3 matrix
#'   ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#' 
#'   # Create an initial stage structure
#'   ( initial <- c(1,3,2))
#' 
#'   # Calculate the matrix of transfer functions for the upper bound on
#'   # inertia, using default arguments
#'   ( tfmat<-tfam_inertia(A,bound="upper") )
#' 
#'   # Plot the transfer function using the S3 method (defaults to p 
#'   # and inertia in this case)
#'   plot(tfmat)
#' 
#'   # Plot inertia against lambda using the S3 method
#'   plot(tfmat, xvar="lambda", yvar="inertia")
#'     
#'   # Plot the transfer function of element [3,2] without the S3 method
#'   par(mfrow=c(1,1))
#'   par(mar=c(5,4,4,2)+0.1)
#'   plot(tfmat$inertia[3,2,]~tfmat$p[3,2,],xlab="p",ylab="lambda",type="l")
#' 
#'   # Create a new matrix with fission of adults
#'   B <- A; B[2,3] <- 0.9; B
#' 
#'   # Calculate the matrix of transfer functions for specified
#'   # initial stage structure, using chosen arguments
#'   # that give the exact structure of the new matrix
#'   # and perturb a minimum of half the value of an element and
#'   # a maximum of double the value of an element
#'   ( etype <- matrix(c(NA, "F", "F", "P", "P", "F", NA, "P", "P"), 
#'                     ncol=3, byrow=TRUE) )
#'   ( tfmat2 <- tfam_inertia(B, vector=initial, elementtype=etype, 
#'                                 Flim=c(-0.5,2), Plim=c(-0.5,2)) )
#' 
#'   # Plot the new matrix of transfer functions using the S3 method
#'   plot(tfmat2)
#'     
#' @concept 
#' KEYWORDS
#'
#' @export
#'
tfam_inertia<-function(A, bound=NULL, vector="n", elementtype=NULL, Flim=c(-1,10), Plim=c(-1,10), plength=100, digits=1e-10){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!is.matrix_irreducible(A)) stop("Matrix A is reducible")
if(!is.matrix_primitive(A)) warning("Matrix A is imprimitive")
laymat<-numeric(length(A))
dim(laymat)<-dim(A)
laymat[which(A>0)]<-1
laymat<-t(t(laymat)*cumsum(t(laymat)))
if(!is.null(elementtype)){
    type<-elementtype
}
if(is.null(elementtype)){
    type<-A
    type[type==0]<-NA
    type[1,][!is.na(type[1,])]<-"F"
    type[2:order,][!is.na(type[2:order,])&type[2:order,]<=1]<-"P"
    type[2:order,][!is.na(type[2:order,])&!type[2:order,]=="P"]<-"F"
}
p<-numeric(order*order*plength)
dim(p)<-c(order,order,plength)
lambda<-p
inertia<-p
for(i in 1:order){
    for(j in 1:order){
        if(A[i,j]!=0){
            d<-matrix(0,order)
            d[i,1]<-1
            e<-matrix(0,order)
            e[j,1]<-1
            if(type[i,j]=="P"){
                minp<-Plim[1]*A[i,j]
                maxp<-Plim[2]*A[i,j]-A[i,j]
                if(maxp>(1-A[i,j])) maxp<-(1-A[i,j])
                if(minp<(-A[i,j])) minp<-(-A[i,j])
                pert<-seq(minp,maxp,(maxp-minp)/(plength-1))
            }
            if(type[i,j]=="F"){
                minp<-Flim[1]*A[i,j]
                maxp<-Flim[2]*A[i,j]-A[i,j]
                if(minp<(-A[i,j])) minp<-(-A[i,j])
                pert<-seq(minp,maxp,(maxp-minp)/(plength-1))
            }
            transfer<-tfa_inertia(A,d,e,bound=bound,vector=vector,prange=pert,digits=digits)
            p[i,j,]<-c(transfer$p,rep(NA,plength-length(transfer$p)))
            lambda[i,j,]<-c(transfer$lambda,rep(NA,plength-length(transfer$lambda)))
            inertia[i,j,]<-c(transfer$inertia,rep(NA,plength-length(transfer$inertia)))
        }
    }
}
final<-list(p=p,lambda=lambda,inertia=inertia,layout=laymat)
class(final)<-c("tfam","list")
return(final)
}
