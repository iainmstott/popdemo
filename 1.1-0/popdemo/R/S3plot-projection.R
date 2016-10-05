################################################################################
#' Plot population dynamics
#'
#' @description
#' Plot dynamics of a population matrix projection model
#'
#' @param x an object of class 'projection' created using \code{\link{project}}.
#'
#' @param labs logical: if \code{TRUE}, then lines are automatically labelled 
#' when \code{x} contains more than one projection.
#'
#' @param plottype for projections generated from dirichlet draws (see 
#' \code{\link{project}}), \code{plottype} has two options. \code{"lines"}
#' will plot each projection as a separate line. \code{"contours"} will 
#' plot shaded contours showing the probabilities of population 
#' densities over over time, calculated across the set of projections from 
#' dirichlet draws. Darker colours indicate more likely population 
#' densities.
#'
#' @param ybreaks if \code{plottype="contour"}, gives the number of breaks on
#' the y axis for generating the grid for the contour plot. A larger number of
#' breaks means a finer resolution grid for the contour plot.
#'
#' @param contourlevels if \code{plottype="contour"}, gives the number of
#' colour/shading levels to use when generating the contour plot. A larger
#' number of levels means a finer resolution colour/shade on the contour
#' plot of population density.
#'
#' @param show.bounds logical: if \code{plottype="contour"}, indicates whether 
#' to plot the bounds on population density: non-zero contours lie within the 
#' bounds.
#'
#' @param ...  arguments to be passed to methods: see \code{\link{par}} and
#' \code{\link{plot}}.
#'
#' @details 
#' Plots population dynamics (time series of density) for objects of class
#' 'projection' created using \code{\link{project}}.  The method is
#' particularly useful for sets of more than one projection.
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
#'   plot(project(A,time=70,standard.A=TRUE,standard.vec=TRUE),
#'        log="y",ylab="Years",xlab="Density")
#' 
#'   # plot a projection of a specified initial stage structure
#'   # for 10 intervals
#'   plot(project(A,vector=initial,time=50))
#' 
#'   # Load the desert Tortoise matrix
#'   data(Tort)
#'
#'   # Project 500 population vectors from a uniform dirichlet 
#'   # distribution, and plot the density of population sizes
#'   # within the bounds of population density
#'   plot(project(Tort, time=30, vector="diri", draws=500, alpha.draws="unif",
#'                standard.A=TRUE),plottype="contour", show.bounds=TRUE)
#'
#' @concept 
#' population projection
#'
#' @export
#'
plot.projection<-function(x,labs=TRUE, plottype = "lines", ybreaks=100, 
                          contourlevels=100, show.bounds=TRUE,...){
if(!is.list(x)){
    if(plottype=="contour") warning('plottype "contour" only valid for projections with vector="diri",\n  defaulting to plottype="lines" instead')
    plottype <- "lines"
    N <- x
    PopulationDensity <- range(N)
}
if (is.list(x)&!any(names(x)=="Bias")){ 
    if(plottype=="contour") warning('plottype "contour" only valid for projections with vector="diri",\n  defaulting to plottype="lines" instead')
    plottype <- "lines"
    N <- x$N
    PopulationDensity <- range(N)
}
if (is.list(x)&any(names(x)=="Bias")){ 
    N <- x$N
    if(show.bounds){ 
        PopulationDensity <- range(x$Bias)
        Bounds <- t(apply(x$Bias,1,range))
    }
    if(!show.bounds) PopulationDensity <- range(N)
}
if (length(dim(N)) == 0) {
    len <- length(N)
    graphics::plot(0:(len - 1), N, type = "l", ...)
}
if (length(dim(N)) == 2 & plottype == "lines") {
    len <- dim(N)[1]
    TimeInterval <- c(0, (len - 1))
    graphics::plot(TimeInterval, PopulationDensity, type = "n", ...)
    options(warn = -1)
    for (i in 1:dim(N)[2]) {
        graphics::lines(0:(len - 1), N[, i], ...)
        if (labs) 
        graphics::text(len - 1, N[len, i], dimnames(N)[[2]][i], adj = c(1, 
             1), ...)
    }
    options(warn = 0)
}
if (length(dim(N)) == 2 & plottype == "contour") {
    len <- dim(N)[1]
    TimeInterval <- c(0, (len - 1))
    ybr<-cbind(seq(PopulationDensity[1]-0.000001,PopulationDensity[2]+0.000001,length.out=ybreaks+1)[-(ybreaks+1)],seq(PopulationDensity[1]-0.000001,PopulationDensity[2]+0.000001,length.out=ybreaks+1)[-1])
    pgrid<-matrix(0,ybreaks,len)
    for(i in 2:(len)){
        pgrid[,i] <- apply(ybr,1,function(x){length(which(N[i,]>=x[1]&N[i,]<x[2]))/dim(N)[2]})
    }
    graphics::plot(TimeInterval, PopulationDensity, type = "n", ...)
    nlevels <- contourlevels
    levels<-seq(min(pgrid),max(pgrid),length.out=nlevels)
    levelcols<-grDevices::colorRampPalette(c("white", "black"))(nlevels)
    graphics::.filled.contour(0:(len-1),apply(ybr,1,mean),t(pgrid),levels=levels,col=levelcols)
    options(warn = -1)
    if(show.bounds){
        graphics::lines(0:(len - 1), Bounds[, 1], col="black", lwd=2)
        graphics::lines(0:(len - 1), Bounds[, 2], col="black", lwd=2)
    }
    options(warn = 0)
}}
