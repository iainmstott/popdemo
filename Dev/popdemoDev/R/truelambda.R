################################################################################
#' Calculate asymptotic growth
#'
#' @description
#' Calculate the true asymptotic growth of a population matrix projection model
#' from the model projection
#'
#' @param A a square, non-negative numeric matrix of any dimension
#' @param vector (optional) a numeric vector or one-column matrix describing 
#' the age/stage distribution used to calculate the projection.
#' @param accuracy the accuracy with which to determine convergence on 
#' asymptotic growth, expressed as a proportion (see details).
#' @param iterations the maximum number of iterations of the model before the 
#' code breaks.  For slowly-converging models and/or high specified convergence
#' accuracy, this may need to be increased.
#'
#' @details 
#' \code{truelambda} works by simulating the given model and manually determining 
#' growth when convergence to the given \code{accuracy} is reached. Convergence 
#' on an asymptotic growth is deemed to have been reached when the growth of the 
#' model stays within the window determined by \code{accuracy} for 10*s 
#' iterations of the model, with s equal to the dimension of \code{A}. For example, 
#' projection of an 8 by 8 matrix with convergence accuracy of 1e-2 is deemed to 
#' have converged on asymptotic growth when 10*8=80 consecutive iterations of the 
#' model have a growth within 1-1e-2=0.99 (i.e. 99\%) and 1+1e-2=1.01 (i.e. 101\%) 
#' of each other.
#' 
#' If \code{vector} is specified, then the asymptotic growth of the projection of 
#' \code{vector} through \code{A} is returned. If \code{vector="n"} then 
#' asymptotic growths of the set of 'stage-biased' vectors are calculated. These 
#' projections are achieved using a set of standard basis vectors equal in number 
#' to the dimension of \code{A}. These have every element equal to 0, except for 
#' a single element equal to 1, i.e. for a matrix of dimension 3, the set of 
#' stage-biased vectors are: \code{c(1,0,0)}, \code{c(0,1,0)} and 
#' \code{c(0,0,1)}.
#' 
#' Asymptotic growth should be equal to the dominant eigenvalue of the matrix. For
#' non-ergodic models this may not be the case: asymptotic growth will depend on
#' the population structure that's projected. \code{truelambda} provides a means
#' to check what the true asymptotic growth of a non-ergodic model is.
#'
#' @return 
#' If \code{vector} is specified, a numeric vector of length 2 giving the range in 
#' which asymptoticgrowth of the model lies.
#' 
#' If \code{vector} is not specified, a 2-column matrix with each row giving the 
#' range in which asymptotic growth lies for its corresponding stage-biased 
#' projection: the number of rows is equal to the dimension of \code{A}; the first 
#' row is the range when projecting [1,0,0,...], the second entry is the range when 
#' projecting [0,1,0,...], etc.
#'
#' @references
#' \itemize{
#'  \item Stott et al. (2010) Methods Ecol. Evol., 1, 242-252.
#' }
#'
#' @family ConvergenceMeasures
#'
#' @examples
#'   # Create a 3x3 irreducible PPM
#'   ( A <- matrix(c(0,0,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#'
#'   # Create an initial stage structure
#'   ( initial <- c(1,3,2) )
#'
#'   # Calculate the true asymptotic growth of the stage-biased
#'   # projections of A
#'   truelambda(A)
#'
#'   # Calculate the true asymptotic growth of the projection of
#'   # A and initial
#'   truelambda(A, vector=initial)
#'
#'   # Create a 3x3 reducible, nonergodic PPM
#'   B<-A; B[3,2] <- 0; B
#'
#'   # Calculate the true asymptotic growth of the 3 stage-biased
#'   # projections of B
#'   truelambda(B)
#'
#' @concept 
#' asymptotic growth project projection ergodic nonergodic Perron Frobenius
#'
#' @export truelambda
#' @importFrom expm "%^%"
#'
truelambda <-
function(A,vector="n",accuracy=1e-7,iterations=1e+5){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
M<-A
eigvals<-eigen(M)$values
lmax<-which.max(Re(eigvals))
lambda<-Re(eigvals[lmax])
A<-M/lambda
acr<-accuracy
if(!all(acr>=0&acr<1)) stop("accuracy must be between 0 and 1")
if(vector[1]=="n"){
    I<-diag(order)
    options(warn=-1)
    Nt<-project(A,time=10*order)
    options(warn=0)
    lambdat<-numeric(10*order)
    lambdat<-rep(list(lambdat),order)
    lastntminus1<-as.list(numeric(order))
    lastnt<-as.list(numeric(order))
    t<-numeric(order)
    truelamb<-matrix(0,ncol=2,nrow=order)
    for(i in 1:order){
        for(j in 1:((10*order)-1)){
            lambdat[[i]][j]<-Nt[j+1,i]/Nt[j,i]
        }
        lastntminus1[[i]]<-(A%^%((10*order)-1))%*%I[,i]
        lastnt[[i]]<-A%*%lastntminus1[[i]]
        if(sum(lastnt[[i]])<.Machine$double.xmin|sum(lastnt[[i]])>.Machine$double.xmax) stop("projection calculation exceeded max or min normalized floating-point")
        t[i]<-1
        while(!all(lambdat[[i]]>=(lambdat[[i]][1]*(1-acr))&lambdat[[i]]<=(lambdat[[i]][1]*(1+acr)))&t[i]<iterations){
            lambdat[[i]][1:((10*order)-1)]<-lambdat[[i]][2:(10*order)]
            lastntminus1[[i]]<-A%*%lastntminus1[[i]]
            lastnt[[i]]<-A%*%lastnt[[i]]
            lambdat[[i]][10*order]<-sum(lastnt[[i]])/sum(lastntminus1[[i]])
            t[i]<-t[i]+1
            if(sum(lastnt[[i]])<.Machine$double.xmin|sum(lastnt[[i]])>.Machine$double.xmax) stop("projection calculation exceeded max or min normalized floating-point")
        }
        if(!t[[i]]<iterations){
             stop("Model is not converging, try decreasing accuracy or increasing iterations\n(some imprimitive matrices may also not converge)")
        }
        truelamb[i,]<-c(lambdat[[i]][1]*lambda*(1-acr),lambdat[[i]][1]*lambda*(1+acr))
    }
}
else{
    n0<-vector
    vector<-n0/sum(n0)
    options(warn=-1)
    Nt<-project(A,vector=vector,standard.vec=TRUE,time=10*order)
    options(warn=0)
    lambdat<-numeric(10*order)
    for(i in 1:(10*order)){
        lambdat[i]<-Nt[i+1]/Nt[i]
    }
    lastntminus1<-(A%^%((10*order)-1))%*%vector
    lastnt<-A%*%lastntminus1
    if(sum(lastnt)<.Machine$double.xmin|sum(lastnt)>.Machine$double.xmax) stop("projection calculation exceeded max or min normalized floating-point")
    t<-1
    while(!all(lambdat>=(lambdat[1]*(1-acr))&lambdat<=(lambdat[1]*(1+acr)))&t<iterations){
        lambdat[1:((10*order)-1)]<-lambdat[2:(10*order)]
        lastntminus1<-A%*%lastntminus1
        lastnt<-A%*%lastnt
        lambdat[10*order]<-sum(lastnt)/sum(lastntminus1)
        t<-t+1
        if(sum(lastnt)<.Machine$double.xmin|sum(lastnt)>.Machine$double.xmax) stop("projection calculation exceeded max or min normalized floating-point")
    }
    if(!t<iterations){
         stop("Model is not converging, try decreasing accuracy or increasing iterations\n(some imprimitive matrices may also not converge)")
    }
    truelamb<-c(lambdat[1]*lambda*(1-acr),lambdat[1]*lambda*(1+acr))
}
return(print(truelamb,-log10(acr)))
}

