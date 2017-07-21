################################################################################
#' Project population dynamics
#'
#' @description
#' Project dynamics of a specified population matrix projection model.
#'
#' @param A a square, non-negative numeric matrix of any dimension.
#'
#' @param vector a numeric vector or matrix describing 
#' the age/stage distribution(s) used to calculate the projection. If a matrix
#' is used, each column describes a single population vector, so the number
#' of rows of the matrix should be equal to the matrix dimension, whilst the 
#' number of columns gives the number of vectors to project. \code{vector} may
#' also take either "n" (default) to calculate the set of stage-biased projections 
#' (see details), or "diri" to project random population vectors drawn from a 
#' dirichlet distribution (see details). 
#'
#' @param time the number of projection intervals.
#'
#' @param standard.A (optional) if \code{TRUE}, scales the PPM by dividing all 
#' elements by the dominant eigenvalue. This standardises asymptotic dynamics: 
#' the dominant eigenvalue of the scaled \code{A} is 1. Useful for assessing 
#' transient dynamics.
#'
#' @param standard.vec (optional) if \code{TRUE}, standardises each \code{vector} 
#' to sum to 1, by dividing each vector by its sum. Useful for assessing projection
#' relative to initial population size.
#'
#' @param return.vec (optional) if \code{TRUE}, returns the time series of 
#' demographic (st)age vectors as well as overall population size.
#' 
#' @param draws if \code{vector="diri"}, the number of population vectors drawn
#' from dirichlet.
#'
#' @param alpha.draws if \code{vector="diri"}, the alpha values passed to 
#' \code{rdirichlet}: used to bias draws towards or away from a certain population
#' structure.
#' 
#'
#' @details 
#' If \code{vector} is specified, \code{project} will calculate population 
#' dynamics through time by projecting this vector / these vectors through 
#' \code{A}. If multiple vectors are specified, a separate population projection
#' is calculated for each.\cr\cr 
#' If \code{vector="n"} (default), \code{project} will automatically project the set of 
#' 'stage-biased' vectors of \code{A}. Effectively, each vector is a population
#' consisting of all individuals in one stage. These projections are achieved using a 
#' set of standard basis vectors equal in number to the dimension of \code{A}.
#' The vectors have every element equal to 0, except for a single element equal to 1,  
#' i.e. for a matrix of dimension 3, the set of stage-biased vectors are: 
#' \code{c(1,0,0)}, \code{c(0,1,0)} and \code{c(0,0,1)}. Stage-biased projections are 
#' useful for seeing how extreme transient dynamics can be.\cr\cr
#' If \code{vector="diri"}, \code{project} draws random population vectors from 
#' the dirichlet distribution. \code{draws} gives the number of population vectors
#' to draw. \code{alpha.draws} gives the parameters for the dirichlet and can be
#' used to bias the draws towards or away from certain population structures.
#' The default is \code{alpha.draws="unif"}, which passes \code{rep(1,dim)} (where
#' dim is the dimension of the matrix), resulting in an equal probability of 
#' any random population vector. Relative values in the vector give the population
#' structure to focus the distribution on, and the absolute value of the vector
#' entries (and their sum) gives the strength of the distribution: values greater
#' than 1 make it more likely to draw from nearby that population structure, 
#' whilst values less than 1 make it less likely to draw from nearby that population
#' structure.\cr\cr
#' Projections returned are of length \code{time+1}, as the first element 
#' represents the population at \code{t=0}.\cr\cr
#' Projections have their own S3 plotting method \code{\link{plot.projection}}
#' to enable easy graphing.
#'
#' @return 
#' If \code{vector} is specified, a numeric vector of population sizes of 
#' length \code{time+1} (if a single vector is given), or a numeric matrix 
#' of population projections where each column represents a single population 
#' projection and is of length \code{time+1} (if multiple vectors are given).\cr\cr
#' If \code{vector="n"}, a numeric matrix of population projections where each column 
#' represents a single stage-biased projection and is of length \code{time+1}.\cr\cr
#' If \code{vector="diri"}, a numeric matrix of population projections where each 
#' column represents projection of a single vector draw and each column is of 
#' length \code{time+1}\cr\cr
#' If \code{return.vec=TRUE}, a list with components:
#' \describe{
#' \item{N}{
#' the numeric vector or matrix of population sizes, as above
#' }
#' \item{vec}{
#' If a single \code{vector} is specified, a numeric matrix of demographic 
#' vectors from projection of \code{vector} through \code{A}. Each column 
#' represents the densities of one life stage in the projection.\cr
#' If multiple \code{vector}s are specified, a three-dimensional array of
#' demographic vectors from projection of the set of initial vectors through
#' \code{A}. The first dimension represents time (and is therefore equal to 
#' \code{time+1}). The second dimension represents the densities of 
#' each stage (and is therefore equal to the dimension of \code{A}). 
#' The third dimension represents each individual projection (and is 
#' therefore equal to the number of initial vectors given).\cr
#' If \code{vector="n"}, a three-dimensional array of demographic vectors from
#' projection of the set of stage-biased vectors through \code{A}. The first 
#' dimension represents time (and is therefore equal to \code{time+1}). The 
#' second dimension represents the densities of each stage (and 
#' is therefore equal to the dimension of \code{A}). The third 
#' dimension represents each individual stage-biased projection (and is 
#' therefore also equal to the dimension of \code{A}). \cr
#' If \code{vector="diri"}, a three-dimensional array of demographic vectors from
#' projection of the dirichlet vector draws projected through \code{A}. The first 
#' dimension represents time (and is therefore equal to \code{time+1}). The second 
#' dimension represents the densities of each stage (and 
#' is therefore equal to the dimension of \code{A}). The third 
#' dimension represents projection of each population draw (and is therefore equal
#' to \code{draws}).
#' }
#' }
#' Some examples for understanding the structure of 3D arrays returned when 
#' \code{return.vec=TRUE}: when projecting a 3 by 3 matrix for >10 time intervals 
#' (see examples), element [11,3,2] represents the density of stage 3 at time 10 
#' for either vector 2 (multiple vectors), stage-bias 2 (\code{vector="n"}) or draw 2 
#' (\code{vector="diri"}); note that because element 1 represents t=0, then t=10 
#' is found at element 11. The vector [,3,2] represents the time series of densities 
#' of stage 3 in the projection of vector 2 / stage-bias 2 / draw 2. The matrix [,,2] 
#' represents the time series of all stages in the projection of vector 2 / stage-bias 
#' 2 / draw 2.\cr\cr
#' Note that the projections inherit the labelling from \code{A} and \code{vector}, if
#' it exists. Both stage and vector names are taken from the COLUMN names of \code{A} 
#' and \code{vector} respectively. These may be useful for selecting from the
#' \code{projection} object, and are passed to \code{\link{plot.projection}} for 
#' labelling graphs.
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
#'   ( pr <- project(A, vector="n", time=70) )
#'   plot(pr)
#'
#'   # Select the projection of stage 2 bias
#'   pr[,2]
#'
#'   # Project stage-biased dynamics of standardised A over 30 
#'   # intervals and return demographic vectors
#'   ( pr2 <- project(A, vector="n", time=30, standard.A=TRUE, return.vec=TRUE) )
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
#'   # Create a matrix with four initial stage structures
#'   ( initial4 <- matrix(runif(3*4), 3, 4))
#'   
#'   # Project the 4 initial stage structures for 30 time intervals,
#'   # standardising the vector sizes
#'   ( pr5 <- project(A, vector = initial4, standard.vec=TRUE, time=30))
#'   plot(pr5)
#'   
#'   # Load the desert Tortoise matrix
#'   data(Tort)
#'
#'   # Project 500 population vectors from a uniform dirichlet 
#'   # distribution, and plot the density of population sizes
#'   # within the bounds of population density
#'   pr6 <- project(Tort, time=30, vector="diri", draws=500, alpha.draws="unif",
#'                  standard.A=TRUE)
#'   plot(pr6, bounds=TRUE, plottype="shady")
#'                  
#
#'   
#' @concept 
#' projection project population
#'
#' @export project
#'
project<-
function (A, vector = "n", time = 100, standard.A = FALSE, standard.vec = FALSE, 
          return.vec = FALSE, draws = 1000, alpha.draws = "unif"){
if (any(!is.matrix(A), (dim(A)[1] != dim(A)[2]))) stop("A must be a square matrix")
order <- dim(A)[1]
if (!isIrreducible(A)) {
    warning("Matrix is reducible")
}
else {
    if (!isPrimitive(A)) {
        warning("Matrix is imprimitive")
    }
}
if (standard.A == TRUE) {
    M <- A
    eigvals <- eigen(M)$values
    lmax <- which.max(Re(eigvals))
    lambda <- Re(eigvals[lmax])
    A <- M/lambda
}
I <- diag(order)
VecBias <- numeric((time + 1) * order * order)
dim(VecBias) <- c(time + 1, order, order)
    PopBias <- numeric((time + 1) * order)
    dim(PopBias) <- c(time + 1, order)
    if(is.null(dimnames(A)[[2]])){
        dimnames(VecBias)[[2]] <- paste(rep("Stage", order), 1:order, sep = "")
        dimnames(VecBias)[[3]] <- paste(rep("Bias", order), 1:order, sep = "")
        dimnames(PopBias)[[2]] <- paste(rep("Bias", order), 1:order, sep = "")
    } else{
        dimnames(VecBias)[[2]] <- dimnames(A)[[2]]
        dimnames(VecBias)[[3]] <- paste(rep("Bias ", order), dimnames(A)[[2]], sep = "")
        dimnames(PopBias)[[2]] <- paste(rep("Bias ", order), dimnames(A)[[2]], sep = "")
    }
for (i in 1:order) {
    VecBias[1, , i] <- I[, i]
    PopBias[1, i] <- 1
    for (j in 1:time) {
        VecBias[j + 1, , i] <- A %*% VecBias[j, , i]
        PopBias[j + 1, i] <- sum(VecBias[j + 1, , i])
    }
}
bounds <- t(apply(PopBias,1,range))
if(vector[1] == "n"){
    attr(PopBias,"class") <- c("projection", "matrix")
    if (return.vec) {
        final <- list(N = PopBias, vec = VecBias)
        attr(final,"bounds") <- bounds
        attr(final,"vec") <- "bias"
        attr(final,"class") <- c("projection", "list")
        return(final)
    } else {
        final <- PopBias
        attr(final,"bounds") <- bounds
        attr(final,"vec") <- "bias"
        return(final)
    }
}
if(vector[1] == "diri") {
    if(alpha.draws[1] == "unif") alpha.draws<-rep(1,order)
    if(length(alpha.draws) != order) stop("length of alpha.draws must be equal to matrix dimension")
    VecDraws <- numeric((time + 1) * order * draws)
    dim(VecDraws) <- c(time + 1, order, draws)
    PopDraws <- numeric((time + 1) * draws)
    dim(PopDraws) <- c(time + 1, draws)
    if(is.null(dimnames(A)[[2]])){
        dimnames(VecDraws)[[2]] <- paste(rep("Stage", order), 1:order, sep = "")
        dimnames(VecDraws)[[3]] <- paste(rep("Draw", draws), 1:draws, sep = "")
        dimnames(PopDraws)[[2]] <- paste(rep("Draw", draws), 1:draws, sep = "")
    } else{
        dimnames(VecDraws)[[2]] <- dimnames(A)[[2]]
        dimnames(VecDraws)[[3]] <- paste(rep("Draw", draws), 1:draws, sep = "")
        dimnames(PopDraws)[[2]] <- paste(rep("Draw", draws), 1:draws, sep = "")
    }
    VecDraws[1, , ] <- t(MCMCpack::rdirichlet(draws,alpha.draws))
    PopDraws[1,] <- 1
    for (i in 1:draws) {
        for (j in 1:time) {
            VecDraws[j + 1, , i] <- A %*% VecDraws[j, ,i]
            PopDraws[j + 1, i] <- sum(VecDraws[j + 1, ,i])
        }
    }
    attr(PopDraws,"class") <- c("projection", "matrix")
    if (return.vec) {
        final <- list(N = PopDraws, vec = VecDraws)
        attr(final,"bounds") <- bounds
        attr(final,"vec") <- "diri"
        attr(final,"class") <- c("projection", "list")
        return(final)
    } else {
        final <- PopDraws
        attr(final,"bounds") <- bounds
        attr(final,"vec") <- "diri"
        return(final)
    }        
}
if(vector[1]!="n"&vector[1]!="diri"&(length(vector)%%order)!=0) stop("vector has the wrong dimension(s)")
if(length(vector)>order){
    nvec <- dim(vector)[2]
    n0 <- vector
    if (standard.vec) vector <- apply(n0,2,function(x){x/sum(x)})
    Vec <- numeric((time + 1) * order * nvec)
    dim(Vec) <- c(time + 1, order, nvec)
    Pop <- numeric((time + 1) * nvec)
    dim(Pop) <- c(time + 1, nvec)
    if(is.null(dimnames(A)[[2]]) & is.null(dimnames(vector)[[2]])){
        dimnames(Vec)[[2]] <- paste(rep("Stage", order), 1:order, sep = "")
        dimnames(Vec)[[3]] <- paste(rep("Vector", nvec), 1:nvec, sep = "")
        dimnames(Pop)[[2]] <- paste(rep("Vector", nvec), 1:nvec, sep = "")
    }
    if(is.null(dimnames(A)[[2]]) & !is.null(dimnames(vector)[[2]])){
        dimnames(Vec)[[2]] <- paste(rep("Stage", order), 1:order, sep = "")
        dimnames(Vec)[[3]] <- dimnames(vector)[[2]]
        dimnames(Pop)[[2]] <- dimnames(vector)[[2]]
    }
    if(!is.null(dimnames(A)[[2]]) & is.null(dimnames(vector)[[2]])){
        dimnames(Vec)[[2]] <- dimnames(A)[[2]]
        dimnames(Vec)[[3]] <- paste(rep("Vector", nvec), 1:nvec, sep = "")
        dimnames(Pop)[[2]] <- paste(rep("Vector", nvec), 1:nvec, sep = "")
    }
    if(!is.null(dimnames(A)[[2]]) & !is.null(dimnames(vector)[[2]])){
        dimnames(Vec)[[2]] <- dimnames(A)[[2]]
        dimnames(Vec)[[3]] <- dimnames(vector)[[2]]
        dimnames(Pop)[[2]] <- dimnames(vector)[[2]]
    }
    Vec[1, , ] <- vector
    Pop[1,] <- colSums(vector)
    for (i in 1:nvec) {
        for (j in 1:time) {
            Vec[j + 1, , i] <- A %*% Vec[j, ,i]
            Pop[j + 1, i] <- sum(Vec[j + 1, ,i])
        }
    }
    class(Pop) <- c("projection", "matrix")
    if (return.vec) {
        final <- list(N = Pop, vec = Vec)
        attr(final,"bounds") <- bounds
        attr(final,"vec") <- "multiple"
        attr(final,"class") <- c("projection", "list")
        return(final)
    } else {
        final <- Pop
        attr(final,"bounds") <- bounds
        attr(final,"vec") <- "multiple"
        attr(final,"class") <- c("projection", "list")
        return(final)
    }        
}
if(length(vector)==order){
    n0 <- vector
    if (standard.vec) vector <- n0/sum(n0)
    Vec <- matrix(0, ncol = order, nrow = time + 1)
    Pop <- numeric(time + 1)
    if(is.null(dimnames(A)[[2]])){
        dimnames(Vec)[[2]] <- paste(rep("Stage", order), 1:order, sep = "")
    }
    if(!is.null(dimnames(A)[[2]])){
        dimnames(Vec)[[2]] <- dimnames(A)[[2]]
    }
    Vec[1, ] <- vector
    Pop[1] <- sum(vector)
    for (i in 1:time) {
        Vec[(i + 1), ] <- A %*% Vec[i, ]
        Pop[i + 1] <- sum(Vec[(i + 1), ])
    }
    class(Pop) <- c("projection", "numeric")
    if (return.vec) {
        final <- list(N = Pop, vec = Vec)
        attr(final,"bounds") <- bounds
        attr(final,"vec") <- "single"
        attr(final,"class") <- c("projection", "list")
        return(final)
    } else {
        final <- Pop
        attr(final,"bounds") <- bounds
        attr(final,"vec") <- "single"
        return(final)
    }
}}
