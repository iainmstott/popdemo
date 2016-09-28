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
#' @concept 
#' population projection
#'
#' @export
#'
plot.projection<-function(x,labs=TRUE,...){
if(is.list(x)) N<-x$N else(N<-x)
if(length(dim(N))==0){
    len<-length(N)
    Time.intervals<-0:(len-1)
    Population.density<-N
    graphics::plot(Time.intervals,Population.density,type="l",...)
}
if(length(dim(N))==2){
    len<-dim(N)[1]
    Time.intervals<-c(0,(len-1))
    Population.density<-c(min(N),max(N))
    graphics::plot(Time.intervals,Population.density,type="n",...)
    for(i in 1:dim(N)[2]){
        options(warn=-1)
        graphics::lines(0:(len-1),N[,i],...)
        if(labs) graphics::text(len-1,N[len,i],paste("Bias",i),adj=c(1,1),...)
        options(warn=0)
    }
}}
