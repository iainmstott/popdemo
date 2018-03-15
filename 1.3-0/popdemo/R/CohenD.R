################################################################################
#' Calculate Cohen's cumulative distance
#'
#' @description
#' Calculate Cohen's cumulative distance metric for a population matrix 
#' projection model.
#'
#' @param A a square, irreducible, non-negative numeric matrix of any dimension.
#'
#' @param vector a numeric vector or one-column matrix describing the age/stage 
#' distribution used to calculate the distance.
#'
#' @details 
#' Calculates the cumulative distance metric as outlined in Cohen (1979). 
#' Will not work for reducible matrices and returns a warning for imprimitive 
#' matrices (although will not function for imprimitive matrices with nonzero 
#' imaginary components in the dominant eigenpair).
#'
#' @return 
#' Cohen's D1.
#'
#' @references
#' Cohen (1979) SIAM J. Appl. Math., 36, 169-175.\cr
#' Stott et al. (2011) Ecol. Lett., 14, 959-970.
#'
#' @family DistanceMeasures
#'
#' @examples
#'   # Create a 3x3 PPM
#'   ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#'
#'   # Create an initial stage structure
#'   ( initial <- c(1,3,2) )
#'
#'   # Calculate Cohen cumulative distance
#'   CohenD(A, vector=initial)
#'
#' @concept 
#' distance vector state space
#'
#' @export CohenD
#'
CohenD <-
function(A,vector){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!isIrreducible(A)) stop("Matrix A is reducible")
if(!isPrimitive(A)) warning("Matrix A is imprimitive")
M<-A
reigs<-eigen(M)
leigs<-eigen(t(M))
lmax<-which.max(Re(reigs$values))
lambda<-reigs$values[lmax]
A<-A/lambda
w<-as.matrix(reigs$vectors[,lmax])
v<-as.matrix(leigs$vectors[,lmax])
if(max(Im(w))>0|max(Im(v))>0) stop("Dominant eigenvectors contain nonzero imaginary components")
w<-abs(Re(w))
v<-abs(Re(v))
w<-w/sum(w)
v<-v/as.vector(t(v)%*%w)
vector<-vector/sum(vector)
I<-diag(order)
wv<-w%*%t(v)
D1v<-(solve(I+wv-A)-wv)%*%vector
D1<-sum(abs(D1v))
return(D1)
}

