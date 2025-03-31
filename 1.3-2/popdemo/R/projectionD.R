################################################################################
#' Calculate projection distance
#'
#' @description
#' Calculate projection distance for a population matrix projection model.
#'
#' @param A a square, irreducible, non-negative numeric matrix of any dimension.
#' @param vector a numeric vector or one-column matrix describing the age/stage 
#' distribution used to calculate the distance.
#'
#' @details 
#' \code{projectionD} (Haridas & Tuljapurkar 2007) is the difference 
#' between the reproductive value of a population with demographic distribution 
#' given by \code{vector} and the reproductive value of a population in stable 
#' state.
#' 
#' \code{projectionD} will not work for reducible matrices and returns a 
#' warning for imprimitive matrices (although will not function for imprimitive 
#' matrices with nonzero imaginary components in the dominant eigenpair).
#'
#' @return 
#' Projection distance.
#'
#' @references
#' \itemize{
#'  \item Haridas & Tuljapurkar (2007) Ecol. Lett., 10, 1143-1153.
#'  \item Stott et al. (2011) Ecol. Lett., 14, 959-970.
#' }
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
#'   # Calculate projection distance
#'   projectionD(A, vector=initial)
#'
#' @concept distance vector
#' 
#' @export projectionD
#'
projectionD <-
function(A,vector){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
if(!isIrreducible(A)) stop("Matrix A is reducible")
if(!isPrimitive(A)) warning("Matrix is imprimitive")
M<-A
reigs<-eigen(M)
leigs<-eigen(t(M))
lmax<-which.max(Re(reigs$values))
lambda<-reigs$values[lmax]
w<-as.matrix(reigs$vectors[,lmax])
v<-as.matrix(leigs$vectors[,lmax])
if(max(Im(w))>0|max(Im(v))>0) stop("Dominant eigenvectors contain nonzero imaginary components")
w<-abs(Re(w))
v<-abs(Re(v))
w<-w/sum(w)
v<-v/as.vector(t(v)%*%w)
vector<-vector/sum(vector)
alpha0<-as.vector((t(v)%*%vector)-1)
return(alpha0)
}

