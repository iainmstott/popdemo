### DISPLAY METHODS
#' Plot methods for 'Projection' objects
#' 
#' This page describes print and plot methods for \code{\link{Projection-class}}. 
#' Example code is below, or worked examples using these methods are 
#' available in the "Deterministic population dynamics" and "Stochastic 
#' population dynamics" vignettes.
#' 
#' \describe{
#'  \item{\code{plot}}{ plot a Projection object}
#' }
#'
#' @examples
#'   ### Desert tortoise matrix
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
#'   Tortshade <- project(Tort, time = 30, vector = "diri", standard.A = TRUE,
#'                draws = 500, alpha.draws = "unif")
#'   plot(Tortshade, plottype = "shady", bounds = TRUE)
#'   
#'   ### STOCHASTIC PROJECTIONS
#'   # load polar bear data
#'   data(Pbear)
#'   
#'   # project over 50 years with uniform matrix selection
#'   Pbearvec <- c(0.106, 0.068, 0.106, 0.461, 0.151, 0.108)
#'   p2 <- project(Pbear, Pbearvec, time = 50, Aseq = "unif")
#' 
#'   # stochastic projection information
#'   Aseq(p2)
#'   projtype(p2)
#'   nmat(p2)
#'   
#'   # plot
#'   plot(p2, log = "y")
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
#'   Tortshade <- project(Tort, time = 30, vector = "diri", standard.A = TRUE,
#'                draws = 500, alpha.draws = "unif")
#'   plot(Tortshade, plottype = "shady", bounds = TRUE)
#'   
#'   ### STOCHASTIC PROJECTIONS
#'   # load polar bear data
#'   data(Pbear)
#'   
#'   # project over 50 years with uniform matrix selection
#'   Pbearvec <- c(0.106, 0.068, 0.106, 0.461, 0.151, 0.108)
#'   p2 <- project(Pbear, Pbearvec, time = 50, Aseq = "unif")
#' 
#'   # stochastic projection information
#'   Aseq(p2)
#'   projtype(p2)
#'   nmat(p2)
#'   
#'   # plot
#'   plot(p2, log = "y")
#'
#' @name Projection-plots
#' 
NULL

### DEFINE PARAMETERS
#' @param x an object of class "Projection" generated using \code{\link{project}}
#' @param y not used
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


# PLOT
#' @rdname Projection-plots
#' @export
setGeneric("plot")
#' @rdname Projection-plots
#' @export
setMethod("plot", signature = c(x = "Projection", y = "missing"),
          function(x, y, bounds=FALSE, bounds.args=NULL, labs=TRUE, 
                   plottype = "lines", ybreaks=20, shadelevels=100, ...){
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
              gargs$type <- "n"
              gargs$x <- TimeInterval
              gargs$y <- PopulationDensity
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
                          wraptext(len - 1, N[len, i], labels = dimnames(N)[[2]][i], ...)
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
                                            levels = levels, col = gargs$col)
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
