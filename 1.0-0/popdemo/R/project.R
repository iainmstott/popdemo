################################################################################
#' Project population dynamics
#'
#' @description
#' Project dynamics of a specified population matrix projection model.
#'
#' @param A a square, non-negative numeric matrix of any dimension.
#'
#' @param vector (optional) a numeric vector or one-column matrix describing 
#' the age/stage distribution used to calculate the projection.
#'
#' @param time the number of projection intervals.
#'
#' @param standard.A (optional) if \code{TRUE}, scales the PPM by dividing all 
#' elements by the dominant eigenvalue. This standardises asymptotic dynamics: 
#' the dominant eigenvalue of the scaled \code{A} is 1. Useful for assessing 
#' transient dynamics.
#'
#' @param standard.vec (optional) if \code{TRUE}, standardises \code{vector} to 
#' sum to 1 by scaling it by \code{sum(vector)}. Useful for assessing projection
#' relative to initial population size.
#'
#' @param return.vec (optional) if \code{TRUE}, returns the time series of 
#' vectors of demographic distribution as well as overall population size.
#'
#' @details 
#' If \code{vector} is specified, \code{project} will calculate population 
#' dynamics by projecting this vector through \code{A}.\cr\cr
#' If \code{vector="n"}, \code{project} will automatically project the set of 
#' 'stage-biased' vectors of \code{A}. These projections are achieved using a 
#' set of standard basis vectors equal in number to the dimension of \code{A}.
#' These have every element equal to 0, except for a single element equal to 1,  
#' i.e. for a matrix of dimension 3, the set of stage-biased vectors are: 
#' \code{c(1,0,0)}, \code{c(0,1,0)} and \code{c(0,0,1)}. This is useful for
#' seeing how extreme transient dynamics can be.\cr\cr
#' Projections returned are of length \code{time+1}, as the first element 
#' represents the population at \code{t=0}.\cr\cr
#' Projections have their own S3 plotting method \code{\link{plot.projection}}
#' to enable easy graphing.
#'
#' @return 
#' If \code{vector} is specified, a numeric vector of population sizes of length 
#' \code{time+1}.\cr
#' If \code{vector="n"}, a numeric matrix of population projections: each column 
#' represents a single stage-biased projection and each column is of length 
#' \code{time+1}.\cr\cr
#' If \code{return.vec=TRUE}, a list with components:
#' \describe{
#' \item{N}{
#' the numeric vector or matrix of population sizes, as above
#' }
#' \item{vec}{
#' If \code{vector} is specified, a numeric matrix of demographic vectors from 
#' projection of \code{vector} through \code{A}.  Each column represents the 
#' densities of one life stage in the projection.\cr
#' If \code{vector="n"}, a three-dimensional array of demographic vectors from
#' projection of the set of stage-biased vectors through \code{A}. The first 
#' dimension represents time (and is therefore equal to \code{time+1} in 
#' length). The second dimension represents the densities of each stage (and 
#' is therefore equal to the dimension of \code{A} in length). The third 
#' dimension represents each individual stage-biased projection (and is 
#' therefore also equal to the dimension of \code{A} in length). For example, 
#' when projecting a 3 by 3 matrix for >10 time intervals (see examples), the 
#' density of stage 3 in bias 2 at time 10 is found at element [11,3,2] 
#' (note that because element 1 represents t=0, then t=10 is found at element 11); 
#' the time series of densities of stage 2 in bias 1 is found using [,2,1]; 
#' the matrix of population vectors for bias 2 would be found using [,,2].
#' }
#' }    
#'
#' @seealso
#' \code{\link{plot.projection}}
#'
#' @examples
#'   # Create a 3x3 PPM
#'   ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#'
#'   # Create an initial stage structure
#'   ( initial <- c(1,3,2) )
#'
#'   # Project stage-biased dynamics of A over 70 intervals
#'   ( pr <- project(A,time=70) )
#'   plot(pr)
#'
#'   # Select the projection of stage 2 bias
#'   pr[,2]
#'
#'   # Project stage-biased dynamics of standardised A over 30 
#'   # intervals and return demographic vectors
#'   ( pr2 <- project(A, time=30, standard.A=TRUE, return.vec=TRUE) )
#'   plot(pr2)
#'
#'   #Select the projection of stage 2 bias
#'   pr2$N[,2]
#'
#'   # Select the density of stage 3 in bias 2 at time 10
#'   pr2$vec[11,3,2]
#'
#'   # Select the time series of densities of stage 2 in bias 1
#'   pr2$vec[,2,1]
#'
#'   #Select the matrix of population vectors for bias 2
#'   pr2$vec[,,2]
#'
#'   # Project A over 50 intervals using a specified population structure
#'   ( pr3 <- project(A, vector=initial, time=50) )
#'   plot(pr3)
#'
#'   # Project standardised dynamics of A over 10 intervals using 
#'   # standardised initial structure and return demographic vectors
#'   ( pr4 <- project(A, vector=initial, time=10, standard.vec=TRUE, 
#'                    standard.A=TRUE, return.vec=TRUE) )
#'   plot(pr4)
#'
#'   # Select the time series for stage 1
#'   pr4$vec[,1]
#'   
#' @concept 
#' projection project population
#'
#' @export
#'
project <-
function(A,vector="n",time=100,standard.A=FALSE,standard.vec=FALSE,return.vec=FALSE){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!isIrreducible(A)){
    warning("Matrix is reducible")
}
else{
    if(!isPrimitive(A)) warning("Matrix is imprimitive")
}
if(standard.A==TRUE){
    M<-A
    eigvals<-eigen(M)$values
    lmax<-which.max(Re(eigvals))
    lambda<-Re(eigvals[lmax])
    A<-M/lambda
}
if(vector[1]=="n"){
    I<-diag(order)
    EmptyVec<-matrix(0,nrow=order,ncol=time+1)
    VecBias<-numeric((time+1)*order*order)
    dim(VecBias)<-c(time+1,order,order)
    dimnames(VecBias)[2]<-list(paste(rep("Stage",order),1:order,sep=""))
    dimnames(VecBias)[3]<-list(paste(rep("Bias",order),1:order,sep=""))
    PopBias<-numeric((time+1)*order)
    dim(PopBias)<-c(time+1,order)
    dimnames(PopBias)[2]<-list(paste(rep("Bias",order),1:order,sep=""))
    for (i in 1:order){
        VecBias[1,,i]<-I[,i]
        PopBias[1,i]<-1
            for (j in 1:time){
            VecBias[j+1,,i]<-A%*%VecBias[j,,i]
            PopBias[j+1,i]<-sum(VecBias[j+1,,i])
        }
    }
    class(PopBias)<-c("projection","matrix")
    if(return.vec){
        final<-list(N=PopBias,vec=VecBias)
        class(final)<-c("projection","list")
        return(final)
    }
    else{
        return(PopBias)
    }
}
else{
    n0<-vector
    if(standard.vec){
        vector<-n0/sum(n0)
    }
    Vec<-matrix(0,ncol=order,nrow=time+1)
    Vec[1,]<-vector
    dimnames(Vec)[2]<-list(paste(rep("Stage",order),1:order,sep=""))
    Pop<-numeric(time+1)
    Pop[1]<-sum(vector)
    for(i in 1:time){
        Vec[(i+1),]<-A%*%Vec[i,]
        Pop[i+1]<-sum(Vec[(i+1),])
    }
    class(Pop)<-c("projection","numeric")
    if(return.vec){
        final<-list(N=Pop,vec=Vec)
        class(final)<-c("projection","list")
        return(final)
    }
    else{
        return(Pop)
    }
}}

