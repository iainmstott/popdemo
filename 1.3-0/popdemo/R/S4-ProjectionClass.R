################################################################################
#' 'Projection' object S4 class
#'
#' Projection objects are created using the \code{\link{project}} function. 
#' Primarily, they contain overall population size over time: they can be 
#' treated as a vector (single population projection) or matrix (multiple 
#' population projections; see information on slot ".Data" below). They also 
#' contain further information on the population projection. These extra pieces 
#' of information are described below in the "Slots" section, and the methods 
#' for accessing them appear below. These are:
#' \describe{
#'  \item{\code{vec}}{ access population vectors}
#'  \item{\code{bounds}}{ access bounds on population dynamics}
#'  \item{\code{mat}}{ access projection matrix/matrices used to create projection(s)}
#'  \item{\code{stochSeq}}{access sequence of environmental values/matrices used 
#'  to create projection(s)}
#'  \item{\code{projtype}}{ find out projection type}
#'  \item{\code{vectype}}{ access type of vector used to initiate population projection(s)}
#' }
#' Other methods for accessing basic information from the projection are:
#' \describe{
#'  \item{\code{nproj}}{ access projection matrix/matrices used to create projection}
#'  \item{\code{nmat}}{ number of projection matrices used to create projection(s)}
#'  \item{\code{ntime}}{ number of time intervals}
#' }
#' Plotting and display methods for 'Projection' objects can be found on the 
#' \code{\link{Projection-plots}} page.
#' 
#' In addition to the examples below, see the "Deterministic population dynamics" 
#' and "Stochastic population dynamics" vignettes for worked examples that use 
#' the 'Projection' objects.
#' 
#' @slot .Data One or more time series of population sizes. 
#' 'Projection' objects inherit from a standard array, and can be treated as 
#' such. Therefore, if \code{vector} is specified, the 'Projection' object will 
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
#'  length \code{time+1}.
#' }
#'
#' @slot vec Age- or stage-based population vectors. \code{vec} 
#' will be:
#' \itemize{
#'  \item If a single \code{vector} is specified, a numeric matrix of demographic 
#'  vectors from projection of \code{vector} through \code{A}. Each column 
#'  represents the densities of one life (st)age in the projection.
#'  \item If multiple \code{vector}s are specified, a three-dimensional array of
#'  demographic vectors from projection of the set of initial vectors through
#'  \code{A}. The first dimension represents time (and is therefore equal to 
#'  \code{time+1}). The second dimension represents the densities of 
#'  each stage (and is therefore equal to the dimension of \code{A}). 
#'  The third dimension represents each individual projection (and is 
#'  therefore equal to the number of initial vectors given).
#'  \item If \code{vector="n"}, a three-dimensional array of demographic vectors from
#'  projection of the set of stage-biased vectors through \code{A}. The first 
#'  dimension represents time (and is therefore equal to \code{time+1}). The 
#'  second dimension represents the densities of each stage (and 
#'  is therefore equal to the dimension of \code{A}). The third 
#'  dimension represents each individual stage-biased projection (and is 
#'  therefore also equal to the dimension of \code{A}).
#'  \item If\code{vector="diri"}, a three-dimensional array of demographic vectors from
#'  projection of the dirichlet vector draws projected through \code{A}. The first 
#'  dimension represents time (and is therefore equal to \code{time+1}). The second 
#'  dimension represents the densities of each stage (and 
#'  is therefore equal to the dimension of \code{A}). The third 
#'  dimension represents projection of each population draw (and is therefore equal
#'  to \code{draws}).
#' }
#' 
#' Some examples for understanding the structure of 3D arrays returned when 
#' \code{return.vec=TRUE}: when projecting a 3 by 3 matrix for >10 time intervals, 
#' element [11,3,2] represents the density of stage 3 at time 10 
#' for either vector 2 (multiple vectors), stage-bias 2 (\code{vector="n"}) or draw 2 
#' (\code{vector="diri"}); note that because element 1 represents t=0, then t=10 
#' is found at element 11. The vector [,3,2] represents the time series of densities 
#' of stage 3 in the projection of vector 2 / stage-bias 2 / draw 2. The matrix [,,2] 
#' represents the time series of all stages in the projection of vector 2 / stage-bias 
#' 2 / draw 2.
#' 
#' 
#' Note that the projections inherit the labelling from \code{A} and \code{vector}, if
#' it exists. Both stage and vector names are taken from the COLUMN names of \code{A} 
#' and \code{vector} respectively. These may be useful for selecting from the
#' \code{projection} object, and are used when labelling plots of Projection 
#' objects containing multiple population projections.
#' 
#' 
#' Set \code{return.vec = FALSE} when calling \code{project} to prevent population
#' vectors from being saved: in this case, \code{vec} is equal to 
#' \code{numeric(0)}. This may be necessary when projecting large numbers of
#' vectors, as is the case when \code{vector = "diri"}.
#'
#' @slot bounds The bounds on population dynamics (only for deterministic 
#' projections). These represent the maximum and minimum population sizes 
#' achieveable at each time interval of the projection. \code{bounds} is a 
#' matrix with 2 columns (lower and upper bounds, in that order), and the 
#' number of rows is equal to \code{time + 1}.
#' 
#' @slot mat The matrix/matrices used in the population projection. In their 
#' raw form \code{mat} is always a three-dimensional array, where the third 
#' dimension is used to index the different matrices. However, by using the 
#' \code{mat()} accessor function below, it is possible to choose different ways
#' of representing the matrices (matrix, list, array).
#'
#' @slot stochSeq A list giving the sequence of environmental values or matrices used in the 
#' projection. This will only have \code{length > 1} when doing stochastic by matrix
#' element projections where vital rates depend on multiple external sources of 
#' information (e.g. temperature AND precipitation, competitor density AND nutrient
#'  availability) or a function call is used to generate the environmental sequence. 
#' \itemize{
#'   \item Deterministic projections - where there is only 1 matrix or 
#' environmental value/set of values this will always be \code{rep(1, time)}.
#'   \item Stochastic projections (with more than 1 matrix or a set of environmental
#' values)
#'    \itemize{ 
#'      \item If \code{stochSeq} is passed as a numeric or
#' character vector then this slot will take that value. 
#'      \item If \code{stochSeq} is passed as a matrix describing
#' a random markov process is passed, the\code{stochSeq} slot will be a single 
#' random chain.
#'      \item If \code{stochSeq} is passed as a function call (e.g. describing a
#' distribution of environmental values), then the both the expression and vector
#' of values used will be returned in the \code{stochSeq} slot.
#'  }
#'}
#' @slot projtype The type of projection. Either "deterministic" (single matrix; 
#' time-invariant), or "stochastic" (multiple matrices; time-varying). If "stochastic",
#' \code{projtype} will also provide what type ("stochastic - density depedent",
#' "stochastic - matrix element", or "stochastic - matrix selection"). 
#'
#' @slot vectype The type of vector passed to \code{project}. May be "single" 
#' (one vector; one population projection), "multiple" (more than one vector; 
#' several population projections), "bias" (stage-biased vectors; 
#' \code{vector = "n"}), or "diri" (vectors drawn from the dirichlet 
#' distribution; \code{vector = "diri"}).
#'
#' @seealso 
#' \code{\link{project}} \code{\link{Projection-plots}}
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
#'   vec(pr4)[,1]
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
#' @name Projection-class
#' 
NULL

### DEFINE CLASS
#' @rdname Projection-class
#' @export
Projection <- setClass("Projection", contains = "array",
                       slots = c(vec = "array", bounds = "matrix",
                                 mat = "list", stochSeq = "integer", 
                                 projtype = "character", vectype = "character") )

### CLASS VALIDITY
validProjection <- function(object){
  #get necessary slots
  N <- object@.Data
  errors <- character()
  if(!(length(dim(N)) %in% c(1, 2))){
    .Datamsg <- ".Data dimension must be 2 or 3"
    c(errors, .Datamsg)
  }
  if(!(length(dim(vec(object))) %in% c(2, 3))){
    vecmsg <- "vec dimension must be 2 or 3"
    c(errors, vecmsg)
  }
  if(length(dim(bounds(object))) != 2){
    boundsmsg <- "bounds dimension must be 2"
    c(errors, boundsmsg)
  }
  if(length(dim(mat(object, return = "array"))) != 3){
    matmsg <- "all list elements in mat must be matrices"
  }
  if(!all(stochSeq(object) %in% 1:nmat(object))){
    stochSeqmsg1 <- "all stochSeq must be integer numbers >0 and <nmat"
    c(errors, stochSeqmsg1)
  }
  if(!(projtype(object) %in% c("deterministic", "stochastic"))){
    projtypemsg <- "projtype must be either \"deterministic\" or \"stochastic\""
    c(errors, projtypemsg)
  }
  if(!(vectype(object) %in% c("single", "bias", "diri", "multiple"))){
    vectypemsg <- "vectype must be either \"single\", \"bias\", \"multiple\" or \"diri\""
    c(errors, vectypemsg)
  }
  ifelse(length(errors) == 0, 
         TRUE, 
         stop(errorConstructor(errors), call. = FALSE))
}
setValidity("Projection", validProjection)

### DEFINE GENERICS AND METHODS
#' @rdname Projection-class
#' @param object an object of class "Projection" generated using \code{\link{project}}

# VEC
#' vec: access population vectors
#' @rdname Projection-class
#' @export
setGeneric("vec", function(object){
  standardGeneric("vec")
})
#' @rdname Projection-class
#' @export
setMethod("vec", signature = (object = "Projection"), 
          function(object){
            return(object@vec)
          }
)

# BOUNDS
#' @rdname Projection-class
#' @export
setGeneric("bounds", 
           function(object){
             standardGeneric("bounds")
           }
)
#' @rdname Projection-class
#' @export
setMethod("bounds", signature = (object = "Projection"), 
          function(object){
            return(object@bounds)
          }
)

# MAT
#' @rdname Projection-class
#' @param ... further arguments (see method, below)
#' @export
setGeneric("mat", 
           function(object, ...){
             standardGeneric("mat")
           }
)
#' @rdname Projection-class
#' @param return either "simple", "list", or "array": used for accessing the 'mat'
#' slot from a Projection object. Note that only list or array can be used for 
#' stochastic projections, which have more than one matrix.
#' @export
setMethod("mat", signature = (object = "Projection"),
          function(object, return = "simple"){
            #get necessary slots
            if(length(object@mat) == 0) return(object@mat)
            if(!any(return %in% c("simple", "list", "array"))){
              stop("return must be \"simple\", \"list\" or \"array\"")
            }
            if(return == "simple"){
              if(dim(object@mat)[3] == 1) return(object@mat[,,1])
              if(dim(object@mat)[3] > 1){
                return(lapply(
                  apply(object@mat, 3, list), function(x){x[[1]]}
                ))
              }
            }
            if(return == "list"){
              return(lapply(
                apply(object@mat, 3, list), function(x){x[[1]]}
              ))
            }
            if(return == "array"){
              return(object@mat)
            }
          }
)

# stochSeq
#' @rdname Projection-class
#' @export
setGeneric("stochSeq", 
           function(object){
             standardGeneric("stochSeq")
           }
)
#' @rdname Projection-class
#' @export
setMethod("stochSeq", signature = (object = "Projection"), 
          function(object){
            return(object@stochSeq)
          }
)

# PROJTYPE
#' @rdname Projection-class
#' @export
setGeneric("projtype", 
           function(object){
             standardGeneric("projtype")
           }
)
#' @rdname Projection-class
#' @export
setMethod("projtype", signature = (object = "Projection"), 
          function(object){
            return(object@projtype)
          }
)

# VECTYPE
#' @rdname Projection-class
#' @export
setGeneric("vectype", 
           function(object){
             standardGeneric("vectype")
           }
)
#' @rdname Projection-class
#' @export
setMethod("vectype", signature = (object = "Projection"), 
          function(object){
            return(object@vectype)
          }
)

# NPROJ
#' @rdname Projection-class
#' @export
setGeneric("nproj", 
           function(object){
             standardGeneric("nproj")
           }
)
#' @rdname Projection-class
#' @export
setMethod("nproj", signature = (object = "Projection"), 
          function(object){
            if(length(object@.Data) == 0) return(0)
            if(length(dim(object@.Data)) == 1) return(1)
            if(length(dim(object@.Data)) > 1) return(dim(object@.Data)[2])
          }
)

# NMAT
#' @rdname Projection-class
#' @export
setGeneric("nmat", 
           function(object){
             standardGeneric("nmat")
           }
)
#' @rdname Projection-class
#' @export
setMethod("nmat", signature = (object = "Projection"), 
          function(object){
            if(length(object@mat) == 0) return(0)
            if(length(dim(object@mat)) == 3) return(dim(object@mat)[3])
          }
)

# NTIME
#' @rdname Projection-class
#' @export
setGeneric("ntime", 
           function(object){
             standardGeneric("ntime")
           }
)
#' @rdname Projection-class
#' @export
setMethod("ntime", signature = (object = "Projection"), 
          function(object){
            if(length(object@.Data) == 0) return(0)
            if(length(object@.Data) > 0) return(dim(object@.Data)[1] - 1)
          }
)


### DISPLAY METHODS
#' Display methods for 'Projection' objects
#' 
#' This page describes print and plot methods for \code{\link{Projection-class}}. 
#' Example code is below, or worked examples using these methods are 
#' available in the "Deterministic population dynamics" and "Stochastic 
#' population dynamics" vignettes.
#' 
#' @seealso 
#' \code{\link{project}} \code{\link{Projection-class}}
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
#'   vec(pr4)[,1]
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
#' @name Projection-plots
#' 
NULL

###DEFINE DISPLAY METHODS
#' @rdname Projection-plots
#' 
#' @param object,x an object of class "Projection" generated using \code{\link{project}}
#' 

# SHOW
#' @rdname Projection-plots
#' @export
#' 
setMethod("show", signature = (object = "Projection"),
          function(object){
            #get necessary slots
            N <- object@.Data
            if(length(N) == 0){
              callNextMethod()
            }
            if(length(N) > 0){
              timec <- as.character(length(stochSeq(object)))
              if(length(vectype(object)) > 0){
                if(vectype(object) == "single"){
                  nc1 <- as.character(1)
                  nc2 <- ""
                  vectypec <- ""
                }
                if(vectype(object) == "multiple"){
                  nc1 <- as.character(dim(N)[2])
                  nc2 <- "s"
                  vectypec <- ""
                }
                if(vectype(object) == "bias"){
                  nc1 <- as.character(dim(N)[2])
                  nc2 <- "s"
                  vectypec <- " (stage-biased initial vectors) "
                }
                if(vectype(object) == "diri"){
                  nc1 <- as.character(dim(N)[2])
                  nc2 <- "s"
                  vectypec <- " (dirichlet initial vectors) "
                }
                cat(paste(nc1, " ", projtype(object), " population projection", nc2,
                          vectypec,
                          " over ", timec, " time intervals.\n\n",
                          sep = ""))
                print(N)
              }
            }
          }
)

# PLOT
#' @param bounds logical: indicates whether to plot the bounds on population density.
#' @param bounds.args A list of graphical parameters for plotting the bounds if
#' \code{bounds=T}. The name of each list element indicates the name of the argument.
#' Could include, e.g. \code{list(lwd=2,lty=3,col="darkred")}.
#' @param labs logical: if \code{TRUE}, the plot includes more than one projection
#' and \code{plottype="lines"}, then lines are automatically labelled according to 
#' the names contained in the 'projection' object.
#' @param plottype for projections generated from dirichlet draws (see 
#' \code{\link{project}}), \code{plottype} has two options. \code{"lines"}
#' will plot each projection as a separate line. \code{"shady"} will 
#' plot shaded contours showing the probabilities of population 
#' densities over over time, calculated across the set of projections from 
#' dirichlet draws. By default this shaded plot is a gradient of black to
#' white (with black representing higher probabilities), but this can be
#' overridden by using the 'col' argument (see examples).
#' @param ybreaks if \code{plottype="shady"}, gives the number of breaks on
#' the y axis for generating the grid for the shade plot. A larger number of
#' breaks means a finer resolution grid for the shading.
#' @param shadelevels if \code{plottype="shady"} and a palette of colours
#' is not specified using 'col', then \code{shadelevels} gives the number of
#' colour/shading levels to use when generating the black and white shade plot. 
#' A larger number of levels means a finer resolution on the shade
#' plot of population density (see examples).
#' @param ...  arguments to be passed to methods: see \code{\link{par}} and
#' \code{\link{plot}}.
#' 
#' @export
#' @rdname Projection-plots
#' 
setMethod("plot", signature = (x = "Projection"),
          function(x, 
                   bounds=FALSE,
                   bounds.args=NULL, 
                   labs=TRUE, 
                   plottype = "lines",
                   ybreaks=20, 
                   shadelevels=100, ...) {
            #get necessary slots
            N <- x@.Data
            
            if(plottype == "shady" & vectype(x) != "diri"){
              warning('plottype "shady" only valid for projections generated with vector="diri",\n  defaulting to plottype="lines" instead')
              plottype <- "lines"
            }
            if(bounds & projtype(x) == "stochastic"){
              warning("bounds are not implemented for stochastic projections yet")
              bounds <- FALSE
            }
            if(length(dim(N)) == 1) N1 <- N[1]
            if(length(dim(N)) == 2) N1 <- N[1,]
            if(!bounds) PopulationDensity <- range(N)
            if(bounds) PopulationDensity <- range(N1) * range(bounds(x))
            len <- dim(N)[1]
            TimeInterval <- c(0, (len - 1))
            gargs <- list(...)
            gargs$type <- "n"; gargs$x <- TimeInterval; gargs$y <- PopulationDensity
            if(!("xlab"%in%names(gargs))) gargs$xlab <- "Time Intervals"
            if(!("ylab"%in%names(gargs))) gargs$ylab <- "Population size / density"
            wraplines <- function(..., log, axes, frame.plot, panel.first, panel.last){
              graphics::lines(...)
            }
            wraptext <- function(..., log, axes, frame.plot, panel.first, panel.last, type){
              graphics::text(...)
            }
            do.call(graphics::plot, gargs)
            if(length(dim(N)) == 1) {
              wraplines(0:(len - 1), N, ...)
            }
            if(length(dim(N)) == 2 & plottype == "lines") {
              for (i in 1:dim(N)[2]) {
                wraplines(0:(len - 1), N[, i], ...)
                if (labs){
                  graphics::par(adj=1)
                  wraptext(len - 1, N[len, i], dimnames(N)[[2]][i], ...)
                  graphics::par(adj=0.5)
                }
              }
            }
            if (length(dim(N)) == 2 & plottype == "shady") {
              ybr <- cbind(seq(PopulationDensity[1] - 0.000001, PopulationDensity[2] + 0.000001,
                               length.out=ybreaks + 1)[-(ybreaks + 1)],
                           seq(PopulationDensity[1] - 0.000001,PopulationDensity[2] + 0.000001,
                               length.out=ybreaks + 1)[-1])
              pgrid  <-  matrix(0, ybreaks, len)
              for(i in 2:(len)){
                pgrid[,i] <- apply(ybr, 1, function(x){length(which(N[i,]>=x[1]&N[i,]<x[2]))/dim(N)[2]})
              }
              if(!("col"%in%names(gargs))) {
                levels <- seq(min(pgrid),max(pgrid),length.out=shadelevels)
                gargs$col <- grDevices::colorRampPalette(c("white", "black"))(shadelevels)
              }
              if("col"%in%names(gargs)) {
                levels <- seq(min(pgrid),max(pgrid),length.out=length(gargs$col))
              }
              graphics::.filled.contour(0:(len-1), apply(ybr,1,mean), t(pgrid), 
                                        levels = levels, col=gargs$col)
            }
            if(bounds){
              if(is.null(bounds.args)) bounds.args <- list(col="black", lwd=2)
              args.lwr <- bounds.args
              args.lwr$x <- matrix(c(0:(len - 1), min(N1) * bounds(x)[, 1]), ncol = 2)
              args.upr <- bounds.args
              args.upr$x <- matrix(c(0:(len - 1), max(N1) * bounds(x)[, 2]), ncol = 2)
              do.call(wraplines, args.upr)
              do.call(wraplines, args.lwr)
            }
          }
)

errorConstructor <- function(errors) {
  
  
  n <- seq(1, length(errors), 1)
  c('The following errors were found:\n',
    paste(
      n,
      '. ',
      errors, '
      \n',
      sep = ""))
}
