################################################################################
#' Project population dynamics
#'
#' @description
#' Project dynamics of a specified population matrix projection model.
#'
#' @param ... Arguments passed to individual methods.
#' @param object Either a single square matrix, a list of square matrices,
#' a \code{CompadreMSM}, or \code{CompadreDDM} object. See
#' \code{\link{CompadreDDM}} and \code{\link{CompadreMSM}} for more info on 
#' constructing the latter. 
#' @param vector (optional) a numeric vector or matrix describing 
#' the age/stage distribution(s) used to calculate the projection. Single
#' population vectors can be given either as a numeric vector or 
#' one-column matrix. Multiple vectors are specified as a matrix, where 
#' each column describes a single population vector. Therefore the number
#' of rows of the matrix should be equal to the matrix dimension, whilst the 
#' number of columns gives the number of vectors to project. \code{vector} may
#' also take either "n" (default) to calculate the set of stage-biased projections 
#' (see details), or "diri" to project random population vectors drawn from a 
#' dirichlet distribution (see details). When specifying a numeric vector/matrix for
#' \code{CompadreDDM} objects, take care that
#' the names supplied here match those in \code{matExprs} slot of the \code{object}.
#' For numeric vectors, names are supplied via \code{my_vector = c(name1 = value1,
#' name2 = value2, etc)} whereas for a numeric matrix, row names should be specified 
#' as \code{dimnames(my_matrix)[[1]] <- c(name1, name2,etc)}. 
#' If names are not supplied, they will generated with the following pattern: \code{V1},
#' \code{V2}, \code{V3}, etc for each stage in the matrix. Therefore, expressions
#' in \code{matExprs} would have to include \code{V1}, etc to function properly.
#' @param time the number of projection intervals.
#' @param standard.A (optional) if \code{TRUE}, scales each matrix in \code{A}
#' by dividing all elements by the dominant eigenvalue. This standardises 
#' asymptotic dynamics: the dominant eigenvalue of the scaled matrix is 1. 
#' Useful for assessing transient dynamics. This is not used for 
#' \code{CompadreDDM} objects.
#' @param standard.vec (optional) if \code{TRUE}, standardises each \code{vector} 
#' to sum to 1, by dividing each vector by its sum. Useful for assessing projection
#' relative to initial population size.
#' @param return.vec (optional) if \code{TRUE}, returns the time series of 
#' demographic (st)age vectors as well as overall population size.
#' @param stochSeq (optional, ignored for \code{CompadreDDM} objects) the sequence of 
#' matrices in a stochastic projection. 
#' \code{stochSeq} may be either:
#' \itemize{
#'  \item "unif" (default), which results in every matrix in \code{A} having an 
#'  equal, random chance of being chosen at each timestep.
#'  \item a square, nonnegative left-stochastic matrix describing a first-order 
#'  Markov chain used to choose the matrices. The transitions are defined COLUMNWISE: 
#'  each column j describes the probability of choosing stage (row) i at time t+1, 
#'  given that stage (column) j was chosen at time t. \code{stochSeq}  should have the 
#'  same dimension as the number of matrices in \code{A}. 
#'  \item an integer vector giving a specific sequence which corresponds to the
#'  matrices in \code{A}.
#'  \item a character vector giving a specific sequence which corresponds to the
#'  names of the matrices in \code{A}.
#' }
#' @param Astart (optional, for stochastic projections with \code{CompadreMSM} 
#' only) in a stochastic projection, the matrix with which to
#' initialise the projection (either numeric, corresponding to the matrices in 
#' \code{object}, or character, corresponding to the names of matrices in \code{object}). 
#' When \code{Astart = NULL} (the default), a random initial matrix is chosen.
#' 
#' @param draws if \code{vector="diri"}, the number of population vectors drawn
#' from dirichlet.
#' @param alpha.draws if \code{vector="diri"}, the alpha values passed to 
#' \code{\link[MCMCpack]{rdirichlet}}: used to bias draws towards or away from a certain population
#' structure.
#' @param PREcheck many functions in \code{popdemo} first check Primitivity, 
#' Reducibility and/or Ergodicity of matrices, with associated warnings and/or 
#' errors if a matrix breaks any assumptions. Set \code{PREcheck=FALSE} if you
#' want to bypass these checks. Only used for \code{CompadreMSM} objects.
#' 
#' @details 
#' If \code{vector} is specified, \code{project} will calculate population 
#' dynamics through time by projecting this vector / these vectors through 
#' \code{A}. If multiple vectors are specified, a separate population projection
#' is calculated for each.
#'  
#' If \code{vector="n"}, \code{project} will automatically project the set of 
#' 'stage-biased' vectors of \code{object}. Effectively, each vector is a population
#' consisting of all individuals in one stage. These projections are achieved using a 
#' set of standard basis vectors equal in number to the dimension of \code{object}.
#' The vectors have every element equal to 0, except for a single element equal to 1,  
#' i.e. for a matrix of dimension 3, the set of stage-biased vectors are: 
#' \code{c(1,0,0)}, \code{c(0,1,0)} and \code{c(0,0,1)}. Stage-biased projections are 
#' useful for seeing how extreme transient dynamics can be. When doing this with
#' \code{CompadreDDM} objects, the stage names will automatically be generated with
#' the pattern of \code{V1}, \code{V2}, \code{V}\emph{n}, where \emph{n} is the number
#' of stages. Thus, the expressions in the \code{matExprs} slot of the \code{CompadreDDM}
#' need to contain these!
#' 
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
#' structure. When doing this with
#' \code{CompadreDDM} objects, the stage names will automatically be generated with
#' the pattern of \code{V1}, \code{V2}, \code{V}\emph{n}, where \emph{n} is the number
#' of stages. Thus, the expressions in the \code{matExprs} slot of the \code{CompadreDDM}
#' need to contain these!
#' 
#' Projections returned are of length \code{time+1}, as the first element 
#' represents the population at \code{t=0}.
#' 
#' Projections have their own plotting method (see \code{\link{Projection-plots}})
#' to facilitate graphing.
#' 
#' In addition to the examples below, see the "Deterministic population dynamics" 
#' and "Stochastic population dynamics" vignettes for worked examples that use 
#' the \code{project} function.
#'
#' @return 
#' A \code{\link{Projection-class}} item. 
#' 'Projection' objects inherit from a standard array, and can be treated as 
#' such. Therefore, if if \code{vector} is specified, the 'Projection' object will 
#' behave as: 
#' \itemize{
#'  \item if a single \code{vector} is given, a numeric vector of population sizes 
#'  of length \code{time+1}
#'  \item if multiple \code{vector}s are given, a numeric matrix of population 
#'  projections where each column represents a single population projection and 
#'  is of length \code{time+1}
#'  \item if \code{vector="n"}, a numeric matrix of population projections where each column 
#'  represents a single stage-biased projection and is of length \code{time+1}.
#'  \item if \code{vector="diri"}, a numeric matrix of population projections where each 
#'  column represents projection of a single vector draw and each column is of 
#'  length \code{time+1}
#' }
#' 
#' See documentation on \code{\link{Projection-class}} objects to understand how 
#' to access other slots (e.g. (st)age vectors through the population projection) 
#' and for S4 methods (e.g. plotting projections).
#' Some examples for understanding the structure of 3D arrays returned when 
#' \code{return.vec=TRUE}: when projecting a 3 by 3 matrix for >10 time intervals 
#' (see examples), element [11,3,2] represents the density of stage 3 at time 10 
#' for either vector 2 (multiple vectors), stage-bias 2 (\code{vector="n"}) or draw 2 
#' (\code{vector="diri"}); note that because element 1 represents t=0, then t=10 
#' is found at element 11. The vector [,3,2] represents the time series of densities 
#' of stage 3 in the projection of vector 2 / stage-bias 2 / draw 2. The matrix [,,2] 
#' represents the time series of all stages in the projection of vector 2 / stage-bias 
#' 2 / draw 2.
#' 
#' Note that the projections inherit the labelling from \code{object} and \code{vector}, if
#' it exists. Both stage and vector names are taken from the COLUMN names of \code{A} 
#' and \code{vector} respectively. These may be useful for selecting from the
#' \code{projection} object, and for labelling graphs when plotting Projection
#' objects.
#'
#' @seealso
#' \code{\link{Projection-class}} \code{\link{Projection-plots}} 
#'
#' @examples
#'   ### USING PROJECTION OBJECTS
#' 
#'   # Create a 3x3 PPM
#'   ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#'
#'   # Project stage-biased dynamics of A over 70 intervals
#'   ( pr <- project(A, vector="n", time=70) )
#'   plot(pr)
#' 
#'   # Access other slots
#'   vec(pr)  #time sequence of population vectors
#'   bounds(pr)  #bounds on population dynamics
#'   mat(pr)  #matrix used to create projection
#'   stochSeq(pr)  #sequence of matrices (more useful for stochastic projections)
#'   projtype(pr)  #type of projection
#'   vectype(pr)  #type of vector(s) initiating projection
#'
#'   # Extra information on the projection
#'   nproj(pr)  #number of projections
#'   nmat(pr)  #number of matrices (more usefulk for stochastic projections)
#'   ntime(pr)  #number of time intervals
#'   
#'   # Select the projection of stage 2 bias
#'   pr[,2]
#'
#'   # Project stage-biased dynamics of standardised A over 30 intervals
#'   ( pr2 <- project(A, vector="n", time=30, standard.A=TRUE) )
#'   plot(pr2)
#'
#'   #Select the projection of stage 2 bias
#'   pr2[,2]
#'
#'   # Select the density of stage 3 in bias 2 at time 10
#'   vec(pr2)[11,3,2]
#'
#'   # Select the time series of densities of stage 2 in bias 1
#'   vec(pr2)[,2,1]
#'
#'   #Select the matrix of population vectors for bias 2
#'   vec(pr2)[,,2]
#'
#'   # Create an initial stage structure
#'   ( initial <- c(1,3,2) )
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
#'   vec(pr4)[, 1, 1]
#'
#'   ### DETERMINISTIC PROJECTIONS
#' 
#'   # Load the desert Tortoise matrix
#'   data(Tort)
#'
#'   # Create an initial stage structure
#'   Tortvec1 <- c(8, 7, 6, 5, 4, 3, 2, 1)
#'   
#'   # Create a projection over 30 time intervals
#'   ( Tortp1 <- project(Tort, vector = Tortvec1, time = 10) )
#'
#'   # plot p1
#'   plot(Tortp1)
#'   plot(Tortp1, bounds = TRUE) #with bounds
#'  
#'   # new display parameters
#'   plot(Tortp1, bounds = TRUE, col = "red", bty = "n", log = "y", 
#'        ylab = "Number of individuals (log scale)",
#'        bounds.args = list(lty = 2, lwd = 2) )
#' 
#'   # multiple vectors
#'   Tortvec2 <- cbind(Tortvec1, c(1, 2, 3, 4, 5, 6, 7, 8))
#'   plot(project(Tort, vector = Tortvec2), log = "y")
#'   plot(project(Tort, vector = Tortvec2), log = "y", labs = FALSE) #no labels
#' 
#'   # dirichlet distribution 
#'   # darker shading indicates more likely population size
#'   plot(project(Tort, time = 30, vector = "diri", standard.A = TRUE,
#'                draws = 500, alpha.draws = "unif"),
#'        plottype = "shady", bounds = TRUE)
#'   
#'   ### STOCHASTIC PROJECTIONS
#'   # load polar bear data
#'   data(Pbear)
#'   
#'   # project over 50 years with uniform matrix selection
#'   Pbearvec <- c(0.106, 0.068, 0.106, 0.461, 0.151, 0.108)
#'   p2 <- project(Pbear, Pbearvec, time = 50, stochSeq = "unif")
#' 
#'   # stochastic projection information
#'   stochSeq(p2)
#'   projtype(p2)
#'   nmat(p2)
#'   
#'   # plot
#'   plot(p2, log = "y")
#'   
#' @concept 
#' projection project population
#' 
#' @import methods
#' @export project
#' 
#' @name project
#' 

setGeneric('project', function(object, ...) {
  standardGeneric('project')
})


#' @rdname project 
#' @export
setMethod('project', 
          signature = 'CompadreMSM',
          function (object, vector = "n", time = 100, 
                    standard.A = FALSE, standard.vec = FALSE, 
                    return.vec = TRUE,
                    stochSeq = "unif", Astart = NULL,
                    draws = 1000, alpha.draws = "unif", PREcheck = TRUE) {
            
            A <- object@A
           
            .projectMSM_impl(A = A, 
                             vector = vector,
                             time = time,
                             standard.A = standard.A,
                             standard.vec = standard.vec,
                             return.vec = return.vec, 
                             stochSeq = stochSeq,
                             Astart = Astart, 
                             draws = draws,
                             alpha.draws = alpha.draws,
                             PREcheck = PREcheck)
            }
)

#' @rdname project
#' @export
setMethod('project',
          signature = 'matrix',
          function (object, vector = "n", time = 100, 
                    standard.A = FALSE, standard.vec = FALSE, 
                    return.vec = TRUE,
                    stochSeq = "unif", Astart = NULL,
                    draws = 1000, alpha.draws = "unif", PREcheck = TRUE) {
            
            A <- list(object)
            
            .projectMSM_impl(A = A, 
                             vector = vector,
                             time = time,
                             standard.A = standard.A,
                             standard.vec = standard.vec,
                             return.vec = return.vec, 
                             stochSeq = stochSeq,
                             Astart = Astart, 
                             draws = draws,
                             alpha.draws = alpha.draws,
                             PREcheck = PREcheck)
          }
)

#' @rdname project
#' @export
setMethod('project',
          signature = 'list',
          function (object, vector = "n", time = 100, 
                    standard.A = FALSE, standard.vec = FALSE, 
                    return.vec = TRUE,
                    stochSeq = "unif", Astart = NULL,
                    draws = 1000, alpha.draws = "unif", PREcheck = TRUE) {

            .projectMSM_impl(A = object, 
                             vector = vector,
                             time = time,
                             standard.A = standard.A,
                             standard.vec = standard.vec,
                             return.vec = return.vec, 
                             stochSeq = stochSeq,
                             Astart = Astart, 
                             draws = draws,
                             alpha.draws = alpha.draws,
                             PREcheck = PREcheck)
          }
)

.projectMSM_impl <- function (A, 
                              vector, 
                              time, 
                              standard.A, 
                              standard.vec, 
                              return.vec,
                              stochSeq, Astart,
                              draws, 
                              alpha.draws, 
                              PREcheck) {
  
  if(PREcheck){
    .checkRedAndPrim(A)
  }
  
  errors <- .aPreChecks(A, character())
  ifelse(
    length(errors) > 0,
    stop(.errorConstructor(errors)),
    ""
  )
  
  
  A <- .a2Array(A)
  
  # extract some info about matrices
  matDim <- dim(A)[1]
  stageNames <- dimnames(A)[[2]]
  if(is.null(stageNames)) stageNames <- paste("S", 
                                              as.character(1:matDim),                                                        sep = "")
  nMat <- dim(A)[3]
  matrixNames <- dimnames(A)[[3]]
  if(is.null(matrixNames)) {
    matrixNames <- paste("mat", 1:nMat, sep = '')
  }
  
  # standardisations
  if (standard.A) {
    # now extended for stochastics as well
    A <- standardA(A)
  }
  
  # initialize matrix sequences
  
  .checkstochSeq(stochSeq, matrixNames)
  
  
  MCstart <- .checkAstart(Astart, stochSeq, matrixNames, nMat)
  
  # If, for some reason, they gave an integer vector and the length
  # doesn't match time, just set it to that.
  if(is.integer(stochSeq) & is.null(dim(stochSeq))){
    if(length(stochSeq) != time) time <- length(stochSeq)
  }
  
  # initialize the matrix sequence vector, stage+population vectors,
  # and Projection object
  if(nMat == 1) { 
    MC <- rep(1L, time)
  } else {
    MC <- .initMc(stochSeq, time, MCstart, matrixNames, nMat)
  }
  
  vecData <- .initPopVec(A = A,
                         vector = vector, 
                         time = time,
                         matDim = matDim,
                         stageNames = stageNames,
                         draws = draws,
                         alpha.draws = alpha.draws,
                         standard.vec = standard.vec)
  
  popVec <- vecData$popVec
  vecType <- vecData$vecType
  projType <- vecData$projType
  
  popSum <- .initPopSum(dim(popVec)[3], time, popVec)
  
  # browser()
  out <- .iterateMsm(A, 
                     time, 
                     popVec, 
                     popSum, 
                     MC,
                     vecType, 
                     projType,
                     return.vec)
  
  return(out)
  
}

#' @rdname project
#' @importFrom rlang env_bind !!!
setMethod('project', 
          signature = 'CompadreDDM',
          function(object,
                   vector = NULL,
                   time = 100,
                   standard.A = NULL,
                   standard.vec = FALSE,
                   return.vec = TRUE,
                   stochSeq = NULL,
                   Astart = NULL,
                   draws = 1000, 
                   alpha.draws = 'unif',
                   PREcheck = FALSE
                   ) {

            dataList <- object@dataList
            matExprs <- object@matExprs
            matDim <- object@matExprs$matDim
            
            # uncomment if you need to debug
            # browser()
            
            if(is.null(dataList$popVec)){
              # No user supplied vector
              stageNames <- paste('V', 1:matDim, sep = "")
              popVecData <- .initDDPopVec(vector = vector, 
                                         time = time,
                                         matDim = matDim,
                                         stageNames = stageNames,
                                         draws = draws,
                                         alpha.draws = alpha.draws,
                                         standard.vec = standard.vec)
            } else {
              # User supplied vector
              stageNames <- .getStageNames(dataList$popVec)
              popVecData <- .initDDPopVec(vector = dataList$popVec,
                                          time = time,
                                          matDim = matDim,
                                          stageNames = stageNames,
                                          draws = draws,
                                          alpha.draws = alpha.draws,
                                          standard.vec = standard.vec)
            }
            
            vecOut <- popVecData$popVec
            projType <- popVecData$projType
            vecType <- popVecData$vecType
            
            nVec <- dim(vecOut)[[3]]
            nStages <- dim(vecOut)[[2]]
            
            evalEnv <- .initMatEnv(matExprs, dataList)
            
            # set up place to store outputs of interest. AList is where
            # the list of A's goes.
            AList <- list()
            popSum <- .initPopSum(nVec, time, vecOut)
            
            projType <- 'stochastic'
            bounds <- matrix(c(NA_real_, NA_real_), nrow = 1)
            AListCounter <- 1
            # iterate!
            for(i in seq_len(nVec)) {
              # bind initial vector and place it in the data list
              evalEnv$popVec <- vecOut[1, , i]
              
              rlang::env_bind(evalEnv,
                              !!! vecOut[1, , i])

              for(j in seq(1, time, 1)) {
               
                # iterate the matrix and generate new population vector
                currentMat <- rlang::eval_tidy(evalEnv$matrixExpr)
                .checkCurrentMat(currentMat)
                
                # if it's a good matrix, store it
                AList[[AListCounter]] <- currentMat
                
                newPopVec <- as.numeric(currentMat %*% evalEnv$popVec)
                
                vecOut[j + 1, , i] <- newPopVec
                popSum[j + 1, i] <- sum(newPopVec)
                # assign stage names to new population vector
                names(newPopVec) <- names(evalEnv$popVec)
                
                # bind new names to evalEnv for next iteration
                rlang::env_bind(evalEnv,
                                !!! newPopVec)
                evalEnv$popVec <- newPopVec
                
                AListCounter <- AListCounter + 1
              }
              # stash outputs in the Projection Object
            }
            
            output <- .updateProjOutput(vecOut, 
                                        popSum, 
                                        bounds,
                                        AList, 
                                        time,
                                        MC = NA_integer_,
                                        vecType = vecType, 
                                        projType = projType,
                                        return.vec = return.vec)
            
            
            return(output)
            
          }
)          
