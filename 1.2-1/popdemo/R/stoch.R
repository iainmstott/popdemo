################################################################################
#' Project population dynamics
#'
#' @description
#' Analyse long-term dynamics of a stochastic population matrix projection model.
#'
#' @param A a list of matrices. \code{stoch} uses \code{\link{project}} to 
#' perform a stochastic' projection where the matrix varies with each timestep. 
#' The sequence of matrices is determined using \code{Aseq}. Matrices must be 
#' square, non-negative and numeric, and all matrices must have the same dimension. 
#' 
#' @param what what should be returned. A character vector with possible entries
#' "lambda" (to calcualate stochastic growth), "var" (to calculate variance in 
#' stochastic growth) and/or "all" (to calculate both).
#'
#' @param vector (optional) a numeric vector describing the age/stage 
#' distribution used to calculate the projection. If \code{vector} is not 
#' specified, a random vector is drawn. Long-term stochastic dynamics should 
#' usually be the same for any vector, although if all the matrices in A are 
#' reducible (see \code{\link{isIrreducible}}), that may not be the case.
#'
#' @param Aseq the sequence of matrices in a stochastic projection. 
#' \code{Aseq} may be either:
#' \itemize{
#'  \item "unif" (default), which results in every matrix in \code{A} having an 
#'  equal, random chance of being chosen at each timestep.
#'  \item a square, nonnegative left-stochastic matrix describing a 
#'  first-order markov chain used to choose the matrices. This should have the 
#'  same dimension as the number of matrices in \code{A}. 
#'  \item a numeric vector giving a specific sequence which corresponds to the
#'  matrices in \code{A}.
#'  \item a character vector giving a specific sequence which corresponds to the
#'  names of the matrices in \code{A}.
#' }
#' 
#' @param iterations the number of projection intervals. The default is 1e+5.
#' 
#' @param discard the number of initial projection intervals to discard, to 
#' discount near-term effects arising from the choice of vector. The default is
#' 1e+3
#'
#' @param PREcheck many functions in \code{popdemo} first check Primitivity, 
#' Reducibility and/or Ergodicity of matrices, with associated warnings and/or 
#' errors if a matrix breaks any assumptions. Set \code{PREcheck=FALSE} if you
#' want to bypass these checks.
#'
#' @details 
#' Calculates stochastic growth and its variance for a given stochastic population 
#' matrix projection model.
#'
#' @return 
#' A numeric vector with two possible elements: "lambda" (the stochastic population
#' growth rate) and "var" (the variance in stochastic population growth rate). Values
#' returned depend on what's passed to \code{what}.
#'
#' @examples
#'   # load the Polar bear data
#'   ( data(Pbear) )
#'
#'   # Find the stochastic growth for a time series with uniform probability of each
#'   # matrix
#'   ( lambda_unif <- stoch(Pbear, what = "lambda", Aseq = "unif") )
#'
#'   # Find the variance in stochastic growth for a time series with uniform 
#'   # probability of each matrix
#'   ( var_unif <- stoch(Pbear, what = "var", Aseq = "unif") )
#'                  
#'   # Find stochastic growth and its variance for a time series with a sequence of
#'   # matrices where "bad" years happen with probability q
#'   q <- 0.5
#'   prob_seq <- c(rep(1-q,3)/3, rep(q,2)/2)
#'   Pbear_seq <- matrix(rep(prob_seq,5), 5, 5)
#'   ( var_unif <- stoch(Pbear, what = "var", Aseq = Pbear_seq) )
#
#'   
#' @concept 
#' stochastic growth variance projection project population
#'
#' @export stoch
#' @importClassesFrom markovchain markovchain
#' @import methods
#'
stoch<-
function (A, what = "all", Aseq = "unif", 
          vector = NULL, iterations = 1e+4, discard = 1e+3, 
          PREcheck = TRUE){
##what should be calculated?
ifelse("lambda" %in% what, growth <- TRUE, growth <- FALSE)
ifelse("var" %in% what, variance <- TRUE, variance <- FALSE)
if("all" %in% what){
    growth <- TRUE
    variance <- TRUE
}
if(!growth & !variance) stop('"what" does not contain the right information')
## check structure of A and rearrange into an array
if(!any((is.list(A) & all(sapply(A, is.matrix))), 
        (is.array(A) & length(dim(A)) == 3)) ){
    stop("A must be a list of matrices")
}
if(is.list(A) & length(A) == 1) stop("A must contain two or more matrices")
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
## do PREcheck
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
## get important info from A
order <- dim(A)[1]
stagenames <- dimnames(A)[[2]]
if(is.null(stagenames)) stagenames <- paste("S", as.character(1:order), sep = "")
nmat <- dim(A)[3]
matrixnames <- dimnames(A)[[3]]
## check vector
if(!is.null(vector) & length(vector) != order) stop("vector has the wrong dimension")
if(is.null(vector)) vector <- stats::runif(order)
vector <- vector / sum(vector)
## generate the Markov Chain for projection
if(!any(Aseq[1] == "unif", 
        is.matrix(Aseq), 
        is.numeric(Aseq) & is.null(dim(Aseq)),
        is.character(Aseq) & is.null(dim(Aseq)))){
    stop('Aseq must take "unif", a numeric matrix, a numeric vector, or a character vector')
} #maybe add support for specifying a markovchain object
# uniform probability or transition matrix
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
    MC <- markovchain::rmarkovchain(iterations, MCo, useRCpp = FALSE)
    MC <- as.numeric(MC)
    names(MC) <- matrixnames[MC]
}
# numeric sequence
if(is.numeric(Aseq) & is.null(dim(Aseq))){
    if(min(Aseq) < 1) stop("Entries in Aseq are not all greater than 0")
    if(max(Aseq) > nmat) stop("One or more entries in Aseq are greater than the number of matrices in A")
    if(!all(Aseq%%1 == 0)) stop("One or more entries in Aseq are not integers")
    if(length(Aseq) != time) time <- length(Aseq)
    MC <- Aseq
    names(MC) <- matrixnames[MC]
}
#character sequence (convert to numeric)
if(Aseq[1] != "unif" & is.character(Aseq) & is.null(dim(Aseq))){
    if(!all(Aseq %in% dimnames(A)[[3]])) stop("Names of Aseq aren't all found in A")
    if(length(Aseq) != time) time <- length(Aseq)
    MC <- match(Aseq, dimnames(A)[[3]])
    names(MC) <- matrixnames[MC]
}
##project the model
pr <- project(A = A, vector = vector, time = iterations, Aseq = MC, 
              PREcheck = FALSE)
##calculate the stuff
#work out per-timestep growth
gr <- pr[(discard:iterations) + 1] / pr[discard:iterations]
#find the per-timestep mean growth (stochastic growth)
if(growth) gr_mean <- mean(gr)
#find the per-timestep variance in growth
if(variance) gr_var <- stats::var(gr)
#decompose the growth
if(growth & variance) final <- data.frame(lambda = gr_mean,
                                          var = gr_var)
if(growth & !variance) final <- gr_mean
if(!growth & variance) final <- gr_var
return(final)
}

# if(decomp){
#     #find the per-timestep growth decomposition
#     gr_a <- apply(A[, , MC], 3, eigs, what = "lambda")
#     gr_t <- gr / gr_a
#     gr_decomp <- cbind(asymptotic = gr_a, transient = gr_t)
#     #find the mean growth components
#     if(growth) gr_decomp_mean <- colMeans(gr_decomp)
#     #find the vcv of the decomposed components
#     if(variance) gr_decomp_var <- var(gr_decomp)
#     if(growth & variance) final <- list(lambda = gr_mean, 
#                                         var = gr_var, 
#                                         lambda_decomp = gr_decomp_mean,
#                                         var_decomp = vcv)
#     if(growth & !variance) final <- gr_mean
#     if(!growth & variance) final <- gr_var
# }
