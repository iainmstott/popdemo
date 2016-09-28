################################################################################
#' Calculate sensitivity of inertia using transfer functions
#'
#' @aliases
#' tfs_inertia
#' tfsm_inertia
#'
#' @description
#' Calculate the sensitivity of population inertia of a population matrix 
#' projection model using differentiation of the transfer function.
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
#' @param startval \code{tfs_inertia} calculates the limit of the derivative 
#' of the transfer function as lambda of the perturbed matrix approaches the 
#' dominant eigenvalue of \code{A} (see details). \code{startval} provides a 
#' starting value for the algorithm: the smaller \code{startval} is, the quicker 
#' the algorithm should converge.
#'
#' @param tolerance the tolerance level for determining convergence (see
#' details).
#'
#' @param return.fit if \code{TRUE} (and only if \code{d} and \code{e} are
#' specified), the lambda and sensitivity values obtained
#' from the convergence algorithm are returned alongside the sensitivity at the
#' limit.
#'
#' @param plot.fit if \code{TRUE} then convergence of the algorithm is plotted
#' as sensitivity~lambda.
#'
#' @details 
#' DETAILS
#'
#' @return 
#' For \code{tfs_inertia}, the sensitivity of inertia (or its bound) to 
#' the specified perturbation structure. If \code{return.fit=TRUE} a list 
#' containing components:
#' \describe{
#' \item{sens}{the sensitivity of inertia (or its bound) to the specified 
#' perturbation structure}
#' \item{lambda.fit}{the lambda values obtained in the fitting process}
#' \item{sens.fit}{the sensitivity values obtained in the fitting process.}\cr
#' For \code{tfsm_inertia}, a matrix containing sensitivity of inertia
#' (or its bound) to each separate element of \code{A}.
#' }
#'
#' @references
#' Stott et al. (2012) Methods Ecol. Evol., 3, 673-684.
#'
#' @family TransferFunctionAnalyses
#' @family PerturbationAnalyses
#'
#' @examples
#'   # Create a 3x3 matrix
#'   ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#' 
#'   # Create an initial stage structure    
#'   ( initial <- c(1,3,2) )
#' 
#'   # Calculate the sensitivity matrix for the upper bound on inertia
#'   tfsm_inertia(A, bound="upper",tolerance=1e-7)
#' 
#'   # Calculate the sensitivity of simultaneous perturbation to 
#'   # A[1,3] and A[2,3] for the lower bound on inertia
#'   tfs_inertia(A, d=c(1,0,0), e=c(0,1,1), bound="lower")
#' 
#'   # Calculate the sensitivity of simultaneous perturbation to 
#'   # A[1,3] and A[2,3] for specified initial stage structure
#'   # and return and plot the fitting process
#'   tfs_inertia(A, d=c(1,0,0), e=c(0,1,1), vector=initial,tolerance=1e-7,
#'               return.fit=TRUE,plot.fit=TRUE)
#' 
#' @concept 
#' KEYWORDS
#'
#' @export
#'
tfs_inertia <-
function(A,d=NULL,e=NULL,vector="n",bound=NULL,startval=0.001,tolerance=1e-10,return.fit=FALSE,plot.fit=FALSE){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!is.matrix_irreducible(A)) stop("Matrix is reducible")
if(!is.matrix_primitive(A)) stop("Matrix is imprimitive")
leigs<-eigen(t(A))
lmax<-which.max(Re(leigs$values))
lambda<-Re(leigs$values[lmax])
v<-as.matrix(abs(Re(leigs$vectors[,lmax])))
if(vector[1]=="n"){
    if(!any(bound=="upper",bound=="lower")) stop('Please specify bound="upper", bound="lower" or specify vector')
    n0<-as.matrix(rep(0,order))
    if(bound=="upper") n0[which.max(v),1]<-1
    if(bound=="lower") n0[which.min(v),1]<-1
}
else{
    if(!is.null(bound)) warning("Specification of vector overrides calculation of bound")
    n0<-as.matrix(vector/sum(vector))
}
if(tolerance>startval) stop("tolerance must be smaller than startval")
ones<-as.matrix(rep(1,order))
if(any(is.null(d),is.null(e))) stop("please specify a perturbation structure using d and e")
d<-as.matrix(d)
e<-as.matrix(e)
startlambda<-lambda+startval
f1.1<-.tf(A,b=n0,c=e,z=startlambda)
f1d.1<--.tf(A,b=n0,c=e,z=startlambda,exp=-2)
f2.1<-.tf(A,b=d,c=ones,z=startlambda)
f2d.1<--.tf(A,b=d,c=ones,z=startlambda,exp=-2)
f3.1<-.tf(A,b=d,c=e,z=startlambda,exp=-2)
f3d.1<--2*.tf(A,b=d,c=e,z=startlambda,exp=-3)
f4.1<-.tf(A,b=d,c=e,z=startlambda)^2
f5.1<-.tf(A,b=d,c=e,z=startlambda,exp=-2)
diff1.1<-((f1d.1*f2.1*f3.1)+(f1.1*f2d.1*f3.1)-(f1.1*f2.1*f3d.1))/f3.1^2
diff2.1<-f4.1/f5.1
S1<-diff1.1*diff2.1
eps<-startval/2
limlambda<-lambda+eps
f1.2<-.tf(A,b=n0,c=e,z=limlambda)
f1d.2<--.tf(A,b=n0,c=e,z=limlambda,exp=-2)
f2.2<-.tf(A,b=d,c=ones,z=limlambda)
f2d.2<--.tf(A,b=d,c=ones,z=limlambda,exp=-2)
f3.2<-.tf(A,b=d,c=e,z=limlambda,exp=-2)
f3d.2<--2*.tf(A,b=d,c=e,z=limlambda,exp=-3)
f4.2<-.tf(A,b=d,c=e,z=limlambda)^2
f5.2<-.tf(A,b=d,c=e,z=limlambda,exp=-2)
diff1.2<-((f1d.2*f2.2*f3.2)+(f1.2*f2d.2*f3.2)-(f1.2*f2.2*f3d.2))/f3.2^2
diff2.2<-f4.2/f5.2
S2<-diff1.2*diff2.2
lambdas<-c(startlambda,limlambda)
sensitivities<-c(S1,S2)
while(abs(S1-S2)>tolerance){
    f1.1<-f1.2
    f1d.1<-f1d.2
    f2.1<-f2.2
    f2d.1<-f2d.2
    f3.1<-f3.2
    f3d.1<-f3d.2
    f4.1<-f4.2
    f5.1<-f5.2
    diff1.1<-diff1.2
    diff2.1<-diff2.2
    S1<-S2
    eps<-eps/2
    limlambda<-lambda+eps
    f1.2<-.tf(A,b=n0,c=e,z=limlambda)
    f1d.2<--.tf(A,b=n0,c=e,z=limlambda,exp=-2)
    f2.2<-.tf(A,b=d,c=ones,z=limlambda)
    f2d.2<--.tf(A,b=d,c=ones,z=limlambda,exp=-2)
    f3.2<-.tf(A,b=d,c=e,z=limlambda,exp=-2)
    f3d.2<--2*.tf(A,b=d,c=e,z=limlambda,exp=-3)
    f4.2<-.tf(A,b=d,c=e,z=limlambda)^2
    f5.2<-.tf(A,b=d,c=e,z=limlambda,exp=-2)
    diff1.2<-((f1d.2*f2.2*f3.2)+(f1.2*f2d.2*f3.2)-(f1.2*f2.2*f3d.2))/f3.2^2
    diff2.2<-f4.2/f5.2
    S2<-diff1.2*diff2.2
    lambdas<-c(lambdas,limlambda)
    sensitivities<-c(sensitivities,S2)
}
if(plot.fit) plot(lambdas,sensitivities,xlab="lambda",ylab="sensitivity")
if(!return.fit) return(S2) else(return(list(sens=S2,lambda.fit=lambdas,sens.fit=sensitivities)))
}

tfsm_inertia <-
function(A,vector="n",bound=NULL,startval=0.001,tolerance=1e-10){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!is.matrix_irreducible(A)) stop("Matrix is reducible")
if(!is.matrix_primitive(A)) stop("Matrix is imprimitive")
leigs<-eigen(t(A))
lmax<-which.max(Re(leigs$values))
lambda<-Re(leigs$values[lmax])
v<-as.matrix(abs(Re(leigs$vectors[,lmax])))
if(vector[1]=="n"){
    if(!any(bound=="upper",bound=="lower")) stop('Please specify bound="upper", bound="lower" or specify vector')
    n0<-as.matrix(rep(0,order))
    if(bound=="upper") n0[which.max(v),1]<-1
    if(bound=="lower") n0[which.min(v),1]<-1
}
else{
    if(!is.null(bound)) warning("Specification of vector overrides calculation of bound")
    n0<-as.matrix(vector/sum(vector))
}
if(tolerance>startval) stop("tolerance must be smaller than startval")
ones<-as.matrix(rep(1,order))
S<-matrix(0,order,order)
for(i in 1:order){
    for(j in 1:order){
        if(A[i,j]!=0){
            d<-matrix(0,order)
            d[i,1]<-1
            e<-matrix(0,order)
            e[j,1]<-1
            startlambda<-lambda+startval
            f1.1<-.tf(A,b=n0,c=e,z=startlambda)
            f1d.1<--.tf(A,b=n0,c=e,z=startlambda,exp=-2)
            f2.1<-.tf(A,b=d,c=ones,z=startlambda)
            f2d.1<--.tf(A,b=d,c=ones,z=startlambda,exp=-2)
            f3.1<-.tf(A,b=d,c=e,z=startlambda,exp=-2)
            f3d.1<--2*.tf(A,b=d,c=e,z=startlambda,exp=-3)
            f4.1<-.tf(A,b=d,c=e,z=startlambda)^2
            f5.1<-.tf(A,b=d,c=e,z=startlambda,exp=-2)
            diff1.1<-((f1d.1*f2.1*f3.1)+(f1.1*f2d.1*f3.1)-(f1.1*f2.1*f3d.1))/f3.1^2
            diff2.1<-f4.1/f5.1
            S1<-diff1.1*diff2.1
            eps<-startval/2
            limlambda<-lambda+eps
            f1.2<-.tf(A,b=n0,c=e,z=limlambda)
            f1d.2<--.tf(A,b=n0,c=e,z=limlambda,exp=-2)
            f2.2<-.tf(A,b=d,c=ones,z=limlambda)
            f2d.2<--.tf(A,b=d,c=ones,z=limlambda,exp=-2)
            f3.2<-.tf(A,b=d,c=e,z=limlambda,exp=-2)
            f3d.2<--2*.tf(A,b=d,c=e,z=limlambda,exp=-3)
            f4.2<-.tf(A,b=d,c=e,z=limlambda)^2
            f5.2<-.tf(A,b=d,c=e,z=limlambda,exp=-2)
            diff1.2<-((f1d.2*f2.2*f3.2)+(f1.2*f2d.2*f3.2)-(f1.2*f2.2*f3d.2))/f3.2^2
            diff2.2<-f4.2/f5.2
            S2<-diff1.2*diff2.2
            while(abs(S1-S2)>tolerance){
                f1.1<-f1.2
                f1d.1<-f1d.2
                f2.1<-f2.2
                f2d.1<-f2d.2
                f3.1<-f3.2
                f3d.1<-f3d.2
                f4.1<-f4.2
                f5.1<-f5.2
                diff1.1<-diff1.2
                diff2.1<-diff2.2
                S1<-S2
                eps<-eps/2
                limlambda<-lambda+eps
                f1.2<-.tf(A,b=n0,c=e,z=limlambda)
                f1d.2<--.tf(A,b=n0,c=e,z=limlambda,exp=-2)
                f2.2<-.tf(A,b=d,c=ones,z=limlambda)
                f2d.2<--.tf(A,b=d,c=ones,z=limlambda,exp=-2)
                f3.2<-.tf(A,b=d,c=e,z=limlambda,exp=-2)
                f3d.2<--2*.tf(A,b=d,c=e,z=limlambda,exp=-3)
                f4.2<-.tf(A,b=d,c=e,z=limlambda)^2
                f5.2<-.tf(A,b=d,c=e,z=limlambda,exp=-2)
                diff1.2<-((f1d.2*f2.2*f3.2)+(f1.2*f2d.2*f3.2)-(f1.2*f2.2*f3d.2))/f3.2^2
                diff2.2<-f4.2/f5.2
                S2<-diff1.2*diff2.2
            }
            S[i,j]<-S2
        }
    }
}
return(S)
}
