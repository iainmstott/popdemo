################################################################################
#' Calculate time to convergence
#'
#' @description
#' Calculate the time to convergence of a population matrix projection model
#' from the model projection
#'
#' @param A a square, non-negative numeric matrix of any dimension, 
#' or a CompadreMat object (see RCompadre package).
#' @param vector (optional) a numeric vector or one-column matrix describing 
#' the age/stage distribution used to calculate the projection.
#' @param accuracy the accuracy with which to determine convergence on 
#' asymptotic growth, expressed as a proportion (see details).
#' @param iterations the maximum number of iterations of the model before the 
#' code breaks.  For slowly-converging models and/or high specified convergence
#' accuracy, this may need to be increased.
#'
#' @details 
#' \code{convt} works by simulating the given model and manually 
#' determining growth when convergence to the given \code{accuracy} is reached. 
#' Convergence on an asymptotic growth is deemed to have been reached when the 
#' growth of the model stays within the window determined by \code{accuracy} for 
#' 10*s iterations of the model, with s equal to the dimension of \code{A}. For 
#' example, projection of an 8 by 8 matrix with convergence accuracy of 1e-2 is 
#' deemed to have converged on asymptotic growth when 10*8=80 consecutive 
#' iterations of the model have a growth within 1-1e-2=0.99 (i.e. 99\%) and 
#' 1+1e-2=1.01 (i.e. 101\%) of each other.
#' 
#' If \code{vector} is specified, the convergence time of the projection of 
#' \code{vector} through \code{A} is returned. If \code{vector="n"} then 
#' asymptotic growths of the set of 'stage-biased' vectors are calculated. These 
#' projections are achieved using a set of standard basis vectors equal in number 
#' to the dimension of \code{A}. These have every element equal to 0, except for 
#' a single element equal to 1, i.e. for a matrix of dimension 3, the set of 
#' stage-biased vectors are: \code{c(1,0,0)}, \code{c(0,1,0)} and 
#' \code{c(0,0,1)}.
#' 
#' Due to the way in which convergence is defined, \code{convt} can 
#' only properly work for strongly ergodic models. Therefore, it will not 
#' function for imprimitive (therefore potentially weakly ergodic) or reducible 
#' (therefore potentially nonergodic) models.
#'
#' @return 
#' If \code{vector} is specified, the convergence time of \code{vector} projected
#' through \code{A}.
#' 
#' If \code{vector} is not specified, a numeric vector of convergence times for 
#' corresponding stage-biased projections: the length of the vector is equal to 
#' the dimension of \code{A}; the first entry is the convergence time of 
#' [1,0,0,...], the second entry is the convergence time of [0,1,0,...], etc.).
#'
#' @references
#' \itemize{
#'  \item Stott et al. (2011) Ecol. Lett., 14, 959-970.
#' }
#'
#' @family ConvergenceMeasures
#'
#' @examples
#'  # Create a 3x3 PPM
#'  ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#'
#'  # Create an initial stage structure
#'  ( initial <- c(1,3,2) )
#'
#'  # Calculate the convergence time of the 3 stage-biased 
#'  # populations within 0.1% of lambda-max
#'  ( convt(A, accuracy=1e-3) )
#'
#'  # Calculate the convergence time of the projection of initial and A
#'  # to within 0.001% of lambda-max
#'  ( convt(A, vector=initial, accuracy=1e-5) )
#'
#' @concept 
#' project projection converge convergence
#'
#' @export convt
#' @importFrom expm "%^%"
#' @importClassesFrom RCompadre CompadreMat
#' @importFrom RCompadre matA
#'
convt <-
function(A,vector="n",accuracy=1e-2,iterations=1e+5){
    if(class(A %in% "CompadreMat")){
        A <- matA(A)
    }
    if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
    order<-dim(A)[1]
    if(!isIrreducible(A)) stop("Matrix is reducible")
    if(!isPrimitive(A)) stop("Matrix is imprimitive")
    eigvals<-eigen(A)$values
    lmax<-which.max(Re(eigvals))
    lambda<-Re(eigvals[lmax])
    acr<-accuracy
    if(!all(acr>=0&acr<1)) stop("accuracy must be between 0 and 1")
    if(vector[1]=="n"){
        I<-diag(order)
        Nt<-project(A,time=10*order)
        lambdat<-numeric(10*order)
        lambdat<-rep(list(lambdat),order)
        lastntminus1<-as.list(numeric(order))
        lastnt<-as.list(numeric(order))
        t<-numeric(order)
        for(i in 1:order){
            for(j in 1:(10*order)){
                lambdat[[i]][j]<-Nt[j+1,i]/Nt[j,i]
            }
            lastntminus1[[i]]<-(A%^%((10*order)-1))%*%I[,i]
            lastnt[[i]]<-A%*%lastntminus1[[i]]
            t[i]<-1
            while(!all(lambdat[[i]]>(lambda*(1-acr))&lambdat[[i]]<(lambda*(1+acr)))&t[i]<iterations){
                lambdat[[i]][1:((10*order)-1)]<-lambdat[[i]][2:(10*order)]
                lastntminus1[[i]]<-A%*%lastntminus1[[i]]
                lastnt[[i]]<-A%*%lastnt[[i]]
                lambdat[[i]][10*order]<-sum(lastnt[[i]])/sum(lastntminus1[[i]])
                t[i]<-t[i]+1
                if(sum(lastnt[[i]])<.Machine$double.xmin|sum(lastnt[[i]])>.Machine$double.xmax) stop("Projection calculation exceeded max or min normalized floating-point")
            }
            if(!t[i]<iterations) stop("Model is not converging.  Reduce accuracy or increase iterations")
        }
    }
    else{
        n0<-vector
        vector<-n0/sum(n0)
        Nt<-project(A,vector=vector,time=10*order)
        lambdat<-numeric(10*order)
        for(i in 1:(10*order)){
            lambdat[i]<-Nt[i+1]/Nt[i]
        }
        lastntminus1<-(A%^%((10*order)-1))%*%vector
        lastnt<-A%*%lastntminus1
        t<-1
        while(!all(lambdat>(lambda*(1-acr))&lambdat<(lambda*(1+acr)))&t<iterations){
            lambdat[1:((10*order)-1)]<-lambdat[2:(10*order)]
            lastntminus1<-A%*%lastntminus1
            lastnt<-A%*%lastnt
            lambdat[10*order]<-sum(lastnt)/sum(lastntminus1)
            t<-t+1
            if(sum(lastnt)<.Machine$double.xmin|sum(lastnt)>.Machine$double.xmax) stop("Projection calculation exceeded max or min normalized floating-point")
        }
        if(!t<iterations) stop("Model is not converging.  Reduce accuracy or increase iterations")
    }
return(t)}

