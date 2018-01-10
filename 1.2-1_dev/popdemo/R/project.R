################################################################################
#' Project population dynamics
#'
#' @description
#' Project dynamics of a specified population matrix projection model.
#'
#' @param A a matrix, or list of matrices. If \code{A} is a matrix, then 
#' \code{project} performs a 'deterministic' projection, where the matrix
#' does not change with each timestep. If \code{A} is a list of matrices, then 
#' \code{project} performs a 'stochastic' projection where the matrix varies 
#' with each timestep. The sequence of matrices is determined using \code{Aseq}. 
#' Matrices must be square, non-negative and numeric. If \code{A} is a list, 
#' all matrices must have the same dimension. 'projection' objects inherit
#' names from \code{A}: if \code{A} is a matrix, stage names are inherited from 
#' its column names. If \code{A} is a list, stage names are inherited from the 
#' column names of the first matrix, and matrix names are inherited from the 
#' names of the list elements.
#'
#' @param vector (optional) a numeric vector or matrix describing 
#' the age/stage distribution(s) used to calculate the projection. Single
#' population vectors can be given either as a numeric vector or 
#' one-column matrix. Multiple vectors are specified as a matrix, where 
#' each column describes a single population vector, so the number
#' of rows of the matrix should be equal to the matrix dimension, whilst the 
#' number of columns gives the number of vectors to project. \code{vector} may
#' also take either "n" (default) to calculate the set of stage-biased projections 
#' (see details), or "diri" to project random population vectors drawn from a 
#' dirichlet distribution (see details). 
#'
#' @param time the number of projection intervals.
#'
#' @param standard.A (optional) if \code{TRUE}, scales each matrix in \code{A}
#' by dividing all elements by the dominant eigenvalue. This standardises 
#' asymptotic dynamics: the dominant eigenvalue of the scaled matrix is 1. 
#' Useful for assessing transient dynamics.
#'
#' @param standard.vec (optional) if \code{TRUE}, standardises each \code{vector} 
#' to sum to 1, by dividing each vector by its sum. Useful for assessing projection
#' relative to initial population size.
#'
#' @param return.vec (optional) if \code{TRUE}, returns the time series of 
#' demographic (st)age vectors as well as overall population size.
#' 
#' @param Aseq (optional) the sequence of matrices in a stochastic projection. 
#' \code{Aseq} may be either:
#' \itemize{
#'  \item "unif" (default), which results in every matrix in \code{A} having an 
#'  equal, random chance of being chosen at each timestep.
#'  \item a square, nonnegative left-stochastic matrix describing a 
#'  first-order markov chain used to choose the matrices. This should have the 
#'  same dimension as the number of messages in \code{A}. 
#'  \item a numeric vector giving a specific sequence which corresponds to the
#'  matrices in \code{A}.
#'  \item a character vector giving a specific sequence which corresponds to the
#'  names of the matrices in \code{A}.
#' }
#' 
#' \code{Aseq="unif"} (default), then each matrix in \code{A} has an equal chance 
#' of being chosen at each timestep. If \code{A} is a matrix
#' 
#' @param draws if \code{vector="diri"}, the number of population vectors drawn
#' from dirichlet.
#'
#' @param alpha.draws if \code{vector="diri"}, the alpha values passed to 
#' \code{rdirichlet}: used to bias draws towards or away from a certain population
#' structure.
#' 
#' @param PREcheck many functions in \code{popdemo} first check Primitivity, 
#' Reducibility and/or Ergodicity of matrices, with associated warnings and/or 
#' errors if a matrix breaks any assumptions. Set \code{PREcheck=FALSE} if you
#' want to bypass these checks.
#' 
#'
#' @details 
#' If \code{vector} is specified, \code{project} will calculate population 
#' dynamics through time by projecting this vector / these vectors through 
#' \code{A}. If multiple vectors are specified, a separate population projection
#' is calculated for each.\cr\cr 
#' If \code{vector="n"}, \code{project} will automatically project the set of 
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
#' if\code{vector="diri"}, a three-dimensional array of demographic vectors from
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
#'   # Load the desert Tortoise matrix
#'   data(Tort)
#'
#'   # Project 500 population vectors from a uniform dirichlet 
#'   # distribution, and plot the density of population sizes
#'   # within the bounds of population density
#'   pr5 <- project(Tort, time=30, vector="diri", draws=500, alpha.draws="unif",
#'                  standard.A=TRUE)
#'   plot(pr5)
#'                  
#
#'   
#' @concept 
#' projection project population
#'
#' @export project
#' @importClassesFrom markovchain markovchain
#' @import methods
#'
project<-
function (A, vector = "n", time = 100, 
          standard.A = FALSE, standard.vec = FALSE, return.vec = FALSE,
          Aseq = "unif", 
          draws = 1000, alpha.draws = "unif", PREcheck = TRUE){
if(!any(is.matrix(A), 
        (is.list(A) & all(sapply(A, is.matrix))), 
        (is.array(A) & length(dim(A)) == 3)) ){
    stop("A must be a matrix or list of matrices")
}
if(is.list(A) & length(A) == 1) A <- A[[1]]
if(is.matrix(A)){
    M1 <- A
    dim(A) <- c(dim(A), 1)
    dimnames(A)[[1]] <- dimnames(M1)[[1]]
    dimnames(A)[[2]] <- dimnames(M1)[[2]]
    dimnames(A)[[3]] <- NULL
}
if (is.list(A) & length(A) > 1) {
    numA <- length(A)
    alldim <- sapply(A, dim)
    if (!diff(range(alldim))==0) {
        stop("all matrices in A must be square and have the same dimension as each other")
    }
    dimA <- mean(alldim)
    L <- A
    A <- numeric(dimA * dimA * numA); dim(A) <- c(dimA, dimA, numA)
    for(i in 1:numA){
        A[,,i] <- L[[i]]
    }
    dimnames(A)[[1]] <- dimnames(L[[1]])[[1]]
    dimnames(A)[[2]] <- dimnames(L[[1]])[[2]]
    dimnames(A)[[3]] <- names(L)
}
if(is.array(A) & dim(A)[3] == 1){
    if (dim(A)[1] != dim(A)[2]) stop("A must be a square matrix")
    if(PREcheck){
        if (!isIrreducible(A[,,1])) {
            warning("Matrix is reducible")
        } else {
            if (!isPrimitive(A[,,1])) {
                warning("Matrix is imprimitive")
            }
        }
    }
}
if (is.array(A) & dim(A)[3] > 1) {
    if (dim(A)[1] != dim(A)[2]) stop("all matrices in A must be square")
    if (PREcheck) {
        red <- numeric(0)
        imp <- numeric(0)
        for(i in 1:dim(A)[3]){
            if (!isIrreducible(A[,,i])) {
                red <- c(red, i)
            } else {
                if (!isPrimitive(A[,,i])) {
                    imp <- c(imp, i)
                }
            }
        }
        if (length(red) > 0){
            if(is.null(dimnames(A)[[3]])) red <- as.character(red)
            if(!is.null(dimnames(A)[[3]])) red <- dimnames(A)[[3]][red]
            red <- paste(red, collapse=", ")
            warning(paste(c("One or more matrices are reducible (", red, ")"), 
                          collapse=""))
        }
        if (length(imp) > 0){
            if(is.null(dimnames(A)[[3]])) imp <- as.character(imp)
            if(!is.null(dimnames(A)[[3]])) imp <- dimnames(A)[[3]][imp]
            imp <- paste(imp, collapse=", ")
            warning(paste(c("One or more matrices are imprimitive (", imp, ")"), 
                          collapse=""))
        }
    }
}
order <- dim(A)[1]
stagenames <- dimnames(A)[[2]]
if(is.null(stagenames)) stagenames <- paste("S", as.character(1:order), sep = "")
nmat <- dim(A)[3]
matrixnames <- dimnames(A)[[3]]
if (standard.A == TRUE) {
    MS <- A
    eigvals <- apply(MS, 3, function(x){ eigen(x)$values } )
    lmax <- apply(eigvals, 2, function(x){ which.max(Re(x)) })
    lambda <- numeric(nmat)
    A <- numeric(order * order * nmat); dim(A) <- c(order, order, nmat)
    for(i in 1:nmat){
        lambda[i] <- Re(eigvals[,i][lmax[i]])
        A[,,i] <- MS[,,i]/lambda[i]
    }
}
if (nmat == 1) {
    if(Aseq != "unif") warning("Only 1 matrix in A: ignoring Aseq")
    MC <- rep(1L, time)
    names(MC) <- matrixnames[MC]
}
if (nmat > 1) {
    if(!any(Aseq[1] == "unif", 
            is.matrix(Aseq), 
            is.numeric(Aseq) & is.null(dim(Aseq)),
            is.character(Aseq) & is.null(dim(Aseq)))){
        stop('Aseq must take "unif", a numeric matrix, a numeric vector, or a character vector')
    }
#add support for specifying a markovchain object
    if(any(Aseq[1] == "unif", is.matrix(Aseq))){
        if(Aseq[1] == "unif") MCtm <- matrix(rep(1/nmat, nmat^2), nmat, nmat)
        if(is.matrix(Aseq)) MCtm <- Aseq
        dimnames(MCtm) <- NULL
        if(dim(MCtm)[1] != dim(MCtm)[2]) stop("Aseq is not a square matrix")
        if(dim(MCtm)[1] != nmat){
            stop("Dimensions of Aseq must be equal to number of matrices in A")
        }
        if(!all(colSums(MCtm) == 1)) stop("Columns of Aseq do not sum to 1")
        MCo <- new("markovchain", transitionMatrix = MCtm, byrow = F)
        MC <- markovchain::rmarkovchain(time, MCo, useRCpp = FALSE)
        MC <- as.numeric(MC)
        names(MC) <- matrixnames[MC]
    }
    if(is.numeric(Aseq) & is.null(dim(Aseq))){
        if(min(Aseq) < 1) stop("Entries in Aseq are not all greater than 0")
        if(max(Aseq) > nmat) stop("One or more entries in Aseq are greater than the number of matrices in A")
        if(!all(Aseq%%1 == 0)) stop("One or more entries in Aseq are not integers")
        if(length(Aseq) != time) time <- length(Aseq)
        MC <- Aseq
        names(MC) <- matrixnames[MC]
    }
    if(Aseq[1] != "unif" & is.character(Aseq) & is.null(dim(Aseq))){
        if(!all(Aseq %in% dimnames(A)[[3]])) stop("Names of Aseq aren't all found in A")
        if(length(Aseq) != time) time <- length(Aseq)
        MC <- match(Aseq, dimnames(A)[[3]])
        names(MC) <- matrixnames[MC]
    }
}
I <- diag(order)
VecBias <- numeric((time + 1) * order * order)
dim(VecBias) <- c(time + 1, order, order)
dimnames(VecBias)[[2]] <- stagenames
dimnames(VecBias)[[3]] <- paste("bias", stagenames, sep = "")
PopBias <- numeric((time + 1) * order)
dim(PopBias) <- c(time + 1, order)
dimnames(PopBias)[[2]] <- paste("bias", stagenames, sep = "")
for (i in 1:order) {
    VecBias[1, , i] <- I[, i]
    PopBias[1, i] <- 1
    for (j in 1:time) {
        VecBias[j + 1, , i] <- A[, , MC[j]] %*% VecBias[j, , i]
        PopBias[j + 1, i] <- sum(VecBias[j + 1, , i])
    }
}
bounds <- t(apply(PopBias,1,range))
if(vector[1] == "n"){
    if(nmat == 1){
        attr(PopBias, "bounds") <- bounds
        attr(PopBias, "proj") <- "deterministic"
    }
    if(nmat > 1){
        attr(PopBias, "bounds") <- NA
        attr(PopBias, "proj") <- "stochastic"
    }
    attr(PopBias, "seq") <- MC
    attr(PopBias, "vec") <- "bias"
    attr(PopBias, "class") <- c("projection", "matrix")
    if (return.vec) {
        final <- list(N = PopBias, vec = VecBias)
        if(nmat == 1){
            attr(final, "bounds") <- bounds
            attr(final, "proj") <- "deterministic"
        }
        if(nmat > 1){
            attr(final, "bounds") <- NA
            attr(final, "proj") <- "stochastic"
        }
        attr(final, "seq") <- MC
        attr(final, "vec") <- "bias"
        attr(final, "class") <- c("projection", "list")
        return(final)
    } else {
        final <- PopBias
        return(final)
    }
}
if(vector[1] == "diri") {
    if(alpha.draws[1] == "unif") alpha.draws<-rep(1,order)
    if(length(alpha.draws) != order) stop("length of alpha.draws must be equal to matrix dimension")
    VecDraws <- numeric((time + 1) * order * draws)
    dim(VecDraws) <- c(time + 1, order, draws)
    dimnames(VecDraws)[[2]] <- stagenames
    dimnames(VecDraws)[[3]] <- paste("draw", as.character(1:draws), sep = "")
    PopDraws <- numeric((time + 1) * draws)
    dim(PopDraws) <- c(time + 1, draws)
    dimnames(PopDraws)[[2]] <- paste("draw", as.character(1:draws), sep = "")
    VecDraws[1, , ] <- t(MCMCpack::rdirichlet(draws,alpha.draws))
    PopDraws[1,] <- 1
    for (i in 1:draws) {
        for (j in 1:time) {
            VecDraws[j + 1, , i] <- A[, , MC[j]] %*% VecDraws[j, ,i]
            PopDraws[j + 1, i] <- sum(VecDraws[j + 1, ,i])
        }
    }
    if(nmat == 1){
        attr(PopDraws, "bounds") <- bounds
        attr(PopDraws, "proj") <- "deterministic"
    }
    if(nmat > 1){
        attr(PopDraws, "bounds") <- NA
        attr(PopDraws, "proj") <- "stochastic"
    }
    attr(PopDraws, "seq") <- MC
    attr(PopDraws, "vec") <- "diri"
    attr(PopDraws, "class") <- c("projection", "matrix")
    if (return.vec) {
        final <- list(N = PopDraws, vec = VecDraws)
        if(nmat == 1){
            attr(final, "bounds") <- bounds
            attr(final, "proj") <- "deterministic"
        }
        if(nmat > 1){
            attr(final, "bounds") <- NA
            attr(final, "proj") <- "stochastic"
        }
        attr(final, "seq") <- MC
        attr(final, "vec") <- "diri"
        attr(final, "class") <- c("projection", "list")
        return(final)
    } else {
        final <- PopDraws
        return(final)
    }        
}
if(vector[1]!="n"&vector[1]!="diri"&(length(vector)%%order)!=0) stop("vector has the wrong dimension(s)")
if(length(vector) > order){
    nvec <- dim(vector)[2]
    vectornames <- dimnames(vector)[[2]]
    if(is.null(vectornames)) vectornames <- paste("V", as.character(1:nvec), sep="")
    n0 <- vector
    if (standard.vec) vector <- apply(n0, 2, function(x){x/sum(x)})
    Vec <- numeric((time + 1) * order * nvec)
    dim(Vec) <- c(time + 1, order, nvec)
    dimnames(Vec)[[2]] <- stagenames
    dimnames(Vec)[[3]] <- vectornames
    Pop <- numeric((time + 1) * nvec)
    dim(Pop) <- c(time + 1, nvec)
    dimnames(Pop)[[2]] <- vectornames
    Vec[1, , ] <- vector
    Pop[1,] <- colSums(vector)
    for (i in 1:nvec) {
        for (j in 1:time) {
            Vec[j + 1, , i] <- A[, , MC[j]] %*% Vec[j, ,i]
            Pop[j + 1, i] <- sum(Vec[j + 1, ,i])
        }
    }
    if(nmat == 1){
        attr(Pop, "bounds") <- bounds
        attr(Pop, "proj") <- "deterministic"
    }
    if(nmat > 1){
        attr(Pop, "bounds") <- NA
        attr(Pop, "proj") <- "stochastic"
    }
    attr(Pop, "seq") <- MC
    attr(Pop, "vec") <- "multiple"
    attr(Pop, "class") <- c("projection", "matrix")
    if (return.vec) {
        final <- list(N = Pop, vec = Vec)
        if(nmat == 1){
            attr(final, "bounds") <- bounds
            attr(final, "proj") <- "deterministic"
        }
        if(nmat > 1){
            attr(final, "bounds") <- NA
            attr(final, "proj") <- "stochastic"
        }
        attr(final, "seq") <- MC
        attr(final, "vec") <- "multiple"
        attr(final, "class") <- c("projection", "list")
        return(final)
    } else {
        final <- Pop
        return(final)
    }        
}
if(length(vector)==order){
    n0 <- vector
    if (standard.vec) vector <- n0/sum(n0)
    Vec <- matrix(0, ncol = order, nrow = time + 1)
    Pop <- numeric(time + 1)
    dimnames(Vec)[[2]] <- stagenames
    Vec[1, ] <- vector
    Pop[1] <- sum(vector)
    for (i in 1:time) {
        Vec[(i + 1), ] <- A[, , MC[i]] %*% Vec[i, ]
        Pop[i + 1] <- sum(Vec[(i + 1), ])
    }
    if(nmat == 1){
        attr(Pop, "bounds") <- bounds
        attr(Pop, "proj") <- "deterministic"
    }
    if(nmat > 1){
        attr(Pop, "bounds") <- NA
        attr(Pop, "proj") <- "stochastic"
    }
    attr(Pop, "seq") <- MC
    attr(Pop, "vec") <- "single"
    attr(Pop, "class") <- c("projection", "numeric")
    if (return.vec) {
        final <- list(N = Pop, vec = Vec)
        if(nmat == 1){
            attr(final, "bounds") <- bounds
            attr(final, "proj") <- "deterministic"
        }
        if(nmat > 1){
            attr(final, "bounds") <- NA
            attr(final, "proj") <- "stochastic"
        }
        attr(final, "seq") <- MC
        attr(final, "vec") <- "single"
        attr(final, "class") <- c("projection", "list")
        return(final)
    } else {
        final <- Pop
        return(final)
    }
}}
