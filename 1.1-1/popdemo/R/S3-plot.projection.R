################################################################################
#' Plot population dynamics
#'
#' @description
#' Plot dynamics of a population matrix projection model
#'
#' @param x an object of class 'projection' created using \code{\link{project}}.
#'
#' @param bounds logical: indicates whether to plot the bounds on population density.
#'
#' @param bounds.args A list of graphical parameters for plotting the bounds if
#' \code{bounds=T}. The name of each list element indicates the name of the argument.
#' Could include, e.g. \code{list(lwd=2,lty=3,col="darkred")}.
#'
#' @param labs logical: if \code{TRUE}, the plot includes more than one projection
#' and \code{plottype="lines"}, then lines are automatically labelled according to 
#' the names contained in the 'projection' object.
#'
#' @param plottype for projections generated from dirichlet draws (see 
#' \code{\link{project}}), \code{plottype} has two options. \code{"lines"}
#' will plot each projection as a separate line. \code{"shady"} will 
#' plot shaded contours showing the probabilities of population 
#' densities over over time, calculated across the set of projections from 
#' dirichlet draws. By default this shaded plot is a gradient of black to
#' white (with black representing higher probabilities), but this can be
#' overridden by using the 'col' argument (see examples).
#'
#' @param ybreaks if \code{plottype="shady"}, gives the number of breaks on
#' the y axis for generating the grid for the shade plot. A larger number of
#' breaks means a finer resolution grid for the shading.
#'
#' @param shadelevels if \code{plottype="shady"} and a palette of colours
#' is not specified using 'col', then \code{shadelevels} gives the number of
#' colour/shading levels to use when generating the black and white shade plot. 
#' A larger number of levels means a finer resolution on the shade
#' plot of population density (see examples).
#'
#' @param ...  arguments to be passed to methods: see \code{\link{par}} and
#' \code{\link{plot}}.
#'
#' @details 
#' Plots population dynamics (time series of density) for objects of class
#' 'projection' created using \code{\link{project}}.  
#' 
#' @seealso
#' \code{\link{project}}
#'
#' @examples
#'   # Create a 3x3 PPM
#'   ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#' 
#'   # Create an initial stage structure
#'   ( initial <- c(1,3,2) )
#' 
#'   # plot stage-biased dynamics of A over 70 intervals using
#'   # standardised dynamics and log y axis
#'   plot(project(A, vector="n", time=70, standard.A=TRUE, standard.vec=TRUE),
#'        log="y", xlab="Years", ylab="Density")
#' 
#'   # plot a projection of a specified initial stage structure
#'   # for 50 intervals, using points, removing the box
#'   plot(project(A,vector=initial,time=50), type="p", bty="n")
#' 
#'   # Create a matrix with four initial stage structures
#'   ( initial4 <- matrix(runif(3*4), 3, 4))
#'   
#'   # Project the 4 initial stage structures for 30 time intervals,
#'   # standardising the vector sizes, using red lines
#'   plot(project(A, vector = initial4, standard.vec=TRUE, time=30), col="red")
#'   
#'   # Load the desert Tortoise matrix
#'   data(Tort)
#'
#'   # Project 500 population vectors from a uniform dirichlet 
#'   # distribution, and plot the density of population sizes
#'   # within the bounds of population density
#'   plot(project(Tort, time=30, vector="diri", draws=500, alpha.draws="unif",
#'                standard.A=TRUE),plottype="shady", bounds=TRUE)
#'   
#'   # Plot the tortoise dirichlet projection using log axes and in red
#'   plot(project(Tort, time=30, vector="diri", draws=500, alpha.draws="unif",
#'                standard.A=TRUE), plottype="shady", bounds=TRUE, 
#'        log="y", col=colorRampPalette(c("white","red"))(50),
#'        bounds.args=list(col="red", lwd=2))
#'
#' @concept 
#' population projection
#'
#'@export
#'
plot.projection<-function(x, bounds=FALSE, bounds.args=NULL, labs=TRUE, 
                          plottype = "lines", ybreaks=20, shadelevels=100, ...){
    if(plottype=="shady"&attr(x,"vec")!="diri"){
        warning('plottype "shady" only valid for projections generated with vector="diri",\n  defaulting to plottype="lines" instead')
        plottype <- "lines"
    }
    if(!is.list(x)) N <- x
    if(is.list(x)) N <- x$N
    if(length(dim(N)) == 0) N1 <- N[1]
    if(length(dim(N)) == 2) N1 <- N[1,]
    if(!bounds) PopulationDensity <- range(N)
    if(bounds) PopulationDensity <- range(N1) * range(attr(x,"bounds"))
    if(length(dim(N)) == 0) len <- length(N)
    if(length(dim(N)) == 2) len <- dim(N)[1]
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
    if(length(dim(N)) == 0) {
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
        args.lwr$x <- matrix(c(0:(len - 1), min(N1) * attr(x,"bounds")[, 1]), ncol = 2)
        args.upr <- bounds.args
        args.upr$x <- matrix(c(0:(len - 1), max(N1) * attr(x,"bounds")[, 2]), ncol = 2)
        do.call(wraplines, args.upr)
        do.call(wraplines, args.lwr)
    }}
