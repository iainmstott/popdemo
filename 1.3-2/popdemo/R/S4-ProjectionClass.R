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
#'  \item{\code{Aseq}}{ access projection matrix sequence used to create projection(s)}
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
#' @slot Aseq The sequence of matrices used in the projection. For deterministic
#' projections (where there is only 1 matrix) this will always be \code{rep(1, time)}.
#' For stochastic projections (with more than 1 matrix), if \code{Aseq} is given 
#' to \code{project} as a numeric or character vector then this slot will take 
#' that value. If a matrix describing a random markov process is passed, the 
#' \code{Aseq} slot will be a single random chain.
#'
#' @slot projtype The type of projection. Either "deterministic" (single matrix; 
#' time-invariant), or "stochastic" (multiple matrices; time-varying). 
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
#'   Aseq(pr)  #sequence of matrices (more useful for stochastic projections)
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
#' @name Projection-class
#' 
NULL

### DEFINE CLASS
#' @rdname Projection-class
#' @export
Projection <- setClass("Projection", contains = "array",
                       slots = c(vec = "array", bounds = "matrix",
                                 mat = "array", Aseq = "numeric", 
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
    if(!all(Aseq(object) %in% 1:nmat(object))){
        Aseqmsg1 <- "all Aseq must be integer numbers >0 and <nmat"
        c(errors, Aseqmsg1)
    }
    if(!(projtype(object) %in% c("deterministic", "stochastic"))){
        projtypemsg <- "projtype must be either \"deterministic\" or \"stochastic\""
        c(errors, projtypemsg)
    }
    if(!(vectype(object) %in% c("single", "bias", "diri", "multiple"))){
        vectypemsg <- "vectype must be either \"single\", \"bias\", \"multiple\" or \"diri\""
        c(errors, vectypemsg)
    }
    ifelse(length(errors) == 0, TRUE, errors)
}
setValidity("Projection", validProjection)

### DEFINE GENERICS AND METHODS
#' @param object an object of class "Projection" generated using \code{\link{project}}

# VEC
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

# ASEQ
#' @rdname Projection-class
#' @export
setGeneric("Aseq", 
           function(object){
               standardGeneric("Aseq")
           }
)
#' @rdname Projection-class
#' @export
setMethod("Aseq", signature = (object = "Projection"), 
          function(object){
              return(object@Aseq)
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

# SHOW
#' @rdname Projection-class
#' @export
setGeneric("show")

#' @rdname Projection-class
#' @export
setMethod("show", signature = (object = "Projection"),
          function(object){
              #get necessary slots
              N <- object@.Data
              if(length(N) == 0){
                  callNextMethod()
              }
              if(length(N) > 0){
                  timec <- as.character(length(Aseq(object)))
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

