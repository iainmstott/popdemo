################################################################################
#' Calculate Keyfitz's delta
#'
#' @description
#' Calculate Keyfitz's delta for a population matrix projection model.
#'
#' @param A a square, irreducible, non-negative numeric matrix of any dimension.
#' @param vector a numeric vector or one-column matrix describing the age/stage 
#' distribution used to calculate the distance.
#'
#' @details 
#' Keyfitz's delta is the sum of the differences between the stable demographic 
#' vector (the dominant right eigenvector of \code{A}) and the demographic 
#' distribution vector of the population (given by \code{vector}). 
#' \code{KeyfitzD} will not work for reducible matrices and returns a 
#' warning for imprimitive matrices (although will not function for imprimitive 
#' matrices with nonzero imaginary components in the dominant eigenpair).
#'
#' @return 
#' Keyfitz's delta.
#'
#' @references
#' \itemize{
#'  \item Keyfitz (1968) Introduction to the Mathematics of Populations. Addison-Wesley.
#'  \item Stott et al. (2010) Ecol. Lett., 14, 959-970.
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
#'   # Calculate Keyfitz's delta
#'   KeyfitzD(A, vector=initial)
#'
#' @concept distance
#'
#' @export KeyfitzD
#'
KeyfitzD <-
function(A,vector){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
if(!isIrreducible(A)) stop("Matrix A is reducible")
if(!isPrimitive(A)) warning("Matrix A is imprimitive")
reigs<-eigen(A)
lmax<-which.max(Re(reigs$values))
w<-as.matrix(reigs$vectors[,lmax])
if(max(Im(w))>0) stop("Dominant right eigenvector contains nonzero imaginary components")
w<-abs(Re(w))
w<-w/sum(w)
vector<-as.matrix(vector/sum(vector))
delta<-0.5*norm(vector-w)
return(delta)
}

