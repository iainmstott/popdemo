################################################################################
#' Calculate sensitivity using transfer functions
#'
#' @aliases
#' tfs_lambda
#' tfsm_lambda
#'
#' @description
#' Calculate the sensitivity of the dominant eigenvalue of a population matrix 
#' projection model using differentiation of the transfer function.
#'
#' @param A a square, nonnegative numeric matrix of any dimension.
#'
#' @param d,e numeric vectors that determine the perturbation structure (see
#' details).
#'
#' @param startval \code{tfs_lambda} calculates the limit of the derivative of the
#' transfer function as lambda of the perturbed matrix approaches the dominant 
#' eigenvalue of \code{A} (see details). \code{startval} provides a starting 
#' value for the algorithm: the smaller \code{startval} is, the quicker the 
#' algorithm should converge.
#'
#' @param tolerance the tolerance level for determining convergence (see
#' details).
#'
#' @param return.fit if \code{TRUE} the lambda and sensitivity values obtained
#' from the convergence algorithm are returned alongside the sensitivity at the
#' limit.
#'
#' @param plot.fit if \code{TRUE} then convergence of the algorithm is plotted
#' as sensitivity~lambda.
#'
#' @details 
#' \code{tfs_lambda} and \code{tfsm_lambda} differentiate a transfer function to
#' find sensitivity of the dominant eigenvalue of \code{A} to perturbations. 
#' This provides an alternative method to using matrix eigenvectors to 
#' calculate the sensitivity matrix and is useful as it may incorporate a 
#' greater diversity of perturbation structures.\cr\cr
#' \code{tfs_lambda} evaluates the transfer function of a specific perturbation 
#' structure. The perturbation structure is determined by \code{d\%*\%t(e)}. 
#' Therefore, the rows to be perturbed are determined by \code{d} and the 
#' columns to be perturbed are determined by \code{e}. The values in d and e 
#' determine the relative perturbation magnitude. For example, if only entry
#' [3,2] of a 3 by 3 matrix is to be perturbed, then \code{d = c(0,0,1)} and 
#' \code{e = c(0,1,0)}. If entries [3,2] and [3,3] are to be perturbed with the 
#' magnitude of perturbation to [3,2] half that of [3,3] then \code{d = c(0,0,1)} 
#' and \code{e = c(0,0.5,1)}. \code{d} and \code{e} may also be expressed as 
#' numeric one-column matrices, e.g. \code{d = matrix(c(0,0,1), ncol=1)}, 
#' \code{e = matrix(c(0,0.5,1), ncol=1)}. See Hodgson et al. (2006) for more 
#' information on perturbation structures.\cr\cr
#' \code{tfsm_lambda} returns a matrix of sensitivity values for observed
#' transitions (similar to that obtained when using \code{\link{sens}} to
#' evaluate sensitivity using eigenvectors), where a separate transfer function 
#' for each nonzero element of \code{A} is calculated (each element perturbed 
#' independently of the others).\cr\cr
#' The formula used by \code{tfs_lambda} and \code{tfsm_lambda} cannot be
#' evaluated at lambda-max, therefore it is necessary to find the limit of the
#' formula as lambda approaches lambda-max. This is done using a bisection
#' method, starting at a value of lambda-max + \code{startval}. \code{startval}
#' should be small, to avoid the potential of false convergence. The algorithm
#' continues until successive sensitivity calculations are within an accuracy
#' of one another, determined by \code{tolerance}: a \code{tolerance} of 1e-10
#' means that the sensitivity calculation should be accurate to 10 decimal
#' places. However, as the limit approaches lambda-max, matrices are no longer
#' invertible (singular): if matrices are found to be singular then
#' \code{tolerance} should be relaxed and made larger.\cr\cr
#' For \code{tfs_lambda}, there is an extra option to return and/or plot the above
#' fitting process using \code{return.fit=TRUE} and \code{plot.fit=TRUE}
#' respectively.
#' 
#' @return 
#' For \code{tfs_lambda}, the sensitivity of lambda-max to the specified
#' perturbation structure. If \code{return.fit=TRUE} a list containing
#' components:
#' \describe{
#' \item{sens}{the sensitivity of lambda-max to the specified perturbation 
#' structure}
#' \item{lambda.fit}{the lambda values obtained in the fitting process}
#' \item{sens.fit}{the sensitivity values obtained in the fitting process.}\cr
#' For \code{tfsm_lambda}, a matrix containing sensitivity of lambda-max 
#' to each element of \code{A}.
#' }
#'
#' @references
#' Hodgson et al. (2006) J. Theor. Biol., 70, 214-224.
#'
#' @family TransferFunctionAnalyses
#' @family PerturbationAnalyses
#'
#' @examples
#'   # Create a 3x3 matrix
#'   ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#' 
#'   # Calculate the sensitivity matrix
#'   tfsm_lambda(A)
#' 
#'   # Calculate the sensitivity of simultaneous perturbation to 
#'   # A[1,2] and A[1,3]
#'   tfs_lambda(A, d=c(1,0,0), e=c(0,1,1))
#' 
#'   # Calculate the sensitivity of simultaneous perturbation to 
#'   # A[1,2] and A[1,3] and return and plot the fitting process
#'   tfs_lambda(A, d=c(1,0,0), e=c(0,1,1),
#'              return.fit=TRUE, plot.fit=TRUE)
#' 
#' @concept 
#' KEYWORDS
#'
#' @export tfs_lambda tfsm_lambda
#'
tfs_lambda <-
function(A,d=NULL,e=NULL,startval=0.001,tolerance=1e-10,return.fit=FALSE,plot.fit=FALSE){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!isIrreducible(A)) stop("matrix is reducible")
if(!isPrimitive(A)) warning("matrix is imprimitive")
eigvals<-eigen(A)$values
lmax<-which.max(Re(eigvals))
lambda<-eigvals[lmax]
if(max(Im(lambda))>0) stop("dominant eigenvalue contains nonzero imaginary component")
lambda<-Re(lambda)
if(tolerance>startval) stop("tolerance must be smaller than startval")
if(any(is.null(d),is.null(e))) stop("please specify a perturbation structure using d and e")
d<-as.matrix(d)
e<-as.matrix(e)
startlambda<-lambda+startval
f4.1<-.tf(A,b=d,c=e,z=startlambda)^2
f5.1<-.tf(A,b=d,c=e,z=startlambda,exp=-2)
S1<-f4.1/f5.1
eps<-startval/2
limlambda<-lambda+eps
f4.2<-.tf(A,b=d,c=e,z=limlambda)^2
f5.2<-.tf(A,b=d,c=e,z=limlambda,exp=-2)
S2<-f4.2/f5.2
lambdas<-c(startlambda,limlambda)
sensitivities<-c(S1,S2)
while(abs(S1-S2)>tolerance){
    f4.1<-f4.2
    f5.1<-f5.2
    S1<-S2
    eps<-eps/2
    limlambda<-lambda+eps
    f4.2<-.tf(A,b=d,c=e,z=limlambda)^2
    f5.2<-.tf(A,b=d,c=e,z=limlambda,exp=-2)
    S2<-f4.2/f5.2
    lambdas<-c(lambdas,limlambda)
    sensitivities<-c(sensitivities,S2)
}
if(plot.fit) graphics::plot(lambdas,sensitivities,xlab="lambda",ylab="sensitivity")
if(!return.fit) return(S2) else(return(list(sens=S2,lambda.fit=lambdas,sens.fit=sensitivities)))
}

tfsm_lambda <- 
function(A,startval=0.001,tolerance=1e-10){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!isIrreducible(A)) stop("matrix is reducible")
if(!isPrimitive(A)) warning("matrix is imprimitive")
eigvals<-eigen(A)$values
lmax<-which.max(Re(eigvals))
lambda<-eigvals[lmax]
if(max(Im(lambda))>0) stop("dominant eigenvalue contains nonzero imaginary component")
lambda<-Re(lambda)
if(tolerance>startval) stop("tolerance must be smaller than startval")
S<-matrix(0,order,order)
for(i in 1:order){
    for(j in 1:order){
        if(A[i,j]!=0){
            d<-matrix(0,order)
            d[i,1]<-1
            e<-matrix(0,order)
            e[j,1]<-1
            startlambda<-lambda+startval
            f4.1<-.tf(A,b=d,c=e,z=startlambda)^2
            f5.1<-.tf(A,b=d,c=e,z=startlambda,exp=-2)
            S1<-f4.1/f5.1
            eps<-startval/2
            limlambda<-lambda+eps
            f4.2<-.tf(A,b=d,c=e,z=limlambda)^2
            f5.2<-.tf(A,b=d,c=e,z=limlambda,exp=-2)
            S2<-f4.2/f5.2
            while(abs(S1-S2)>tolerance){
                f4.1<-f4.2
                f5.1<-f5.2
                S1<-S2
                eps<-eps/2
                limlambda<-lambda+eps
                f4.2<-.tf(A,b=d,c=e,z=limlambda)^2
                f5.2<-.tf(A,b=d,c=e,z=limlambda,exp=-2)
                S2<-f4.2/f5.2
            }
        S[i,j]<-S2
        }
    }
}
return(S)
}
