################################################################################
#' Plot transfer function
#'
#' @description
#' Plot a matrix of transfer functions
#'
#' @param x an object of class 'tfam' (transfer function analysis matrix) 
#' created using \code{\link{tfam_lambda}} or \code{\link{tfam_inertia}}.
#'
#' @param xvar,yvar (optional) the variables to plot on the x and y axes. May
#' be \code{"p"}, \code{"lambda"} or \code{"inertia"}. Defaults to
#' \code{xvar="p"} and \code{yvar="lambda"} for objects created using
#' \code{tfam_lambda}, and \code{xvar="p"} and \code{yvar="inertia"} for 
#' objects created using code{tfam_inertia}.
#'
#' @param mar the margin limits on the plots: see \code{\link{par}}
#'
#' @param ...  arguments to be passed to methods: see \code{\link{par}} and
#' \code{\link{plot}}.
#'
#' @details
#' \code{plot.tfam} plots matrices of transfer functions (class 
#' \code{tfam}) created using \code{\link{tfam_lambda}} or 
#' \code{\link{tfam_inertia}}. The plot is laid out to correspond with 
#' the nonzero entries of the matrix used to generate the transfer functions, 
#' for easy visual comparison of how perturbation affects different matrix 
#' elements.
#'
#' @seealso
#' Constructor functions: \code{\link{tfam_lambda}}, \code{\link{tfam_inertia}}
#'
#' @examples
#'   # Create a 3x3 matrix
#'   ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#'
#'   # Calculate the matrix of transfer functions using default arguments
#'   ( tfmat<-tfam_lambda(A) )
#' 
#'   # Plot the matrix of transfer functions
#'   plot(tfmat)
#' 
#'   # Create an initial stage structure
#'   ( initial <- c(1,3,2) )
#'
#'   # Calculate the matrix of transfer functions for inertia and 
#'   # specified initial stage structure using default arguments
#'   ( tfmat2<-tfam_inertia(A,vector=initial) )
#'
#'   # Plot the result (defaults to inertia ~ p)
#'   plot(tfmat2)
#' 
#'   # Plot inertia ~ lambda
#'   plot(tfmat2, xvar="lambda", yvar="inertia")
#' 
#' @concept 
#' transfer function systems control nonlinear perturbation population viability
#' PVA ecology demography PPM MPM
#' 
#' @method plot tfam
#' @export
plot.tfam <- function(x,xvar=NULL,yvar=NULL,mar=c(1.1,1.1,0.1,0.1),...){
    order<-dim(x$p)[1]
    lay<-rbind(rep(max(x$layout)+1,order),x$layout)
    graphics::layout(lay)
    graphics::par(mar=mar)
    if(is.null(xvar)){
        xv<-x$p
        xvar<-"p"
    }
    else{
        if(xvar=="p") xv<-x$p
        if(xvar=="lambda") xv<-x$lambda
        if(xvar=="inertia") xv<-x$inertia
    }
    if(is.null(yvar)){
        yv<-x[[length(x)-1]]
        yvar<-names(x)[length(x)-1]
    }
    else{
        if(yvar=="p") yv<-x$p
        if(yvar=="lambda") yv<-x$lambda
        if(yvar=="inertia") yv<-x$inertia
    }
    elementrow<-x$rows
    elementcol<-x$cols
    for(i in 1:order){
        for(j in 1:order){
            if(x$layout[i,j]!=0){
                graphics::plot(xv[i,j,],yv[i,j,],type="l",xaxt="n",yaxt="n",...)
                graphics::axis(side=1,tck=0.05,padj=-1.7,cex.axis=0.9)
                graphics::axis(side=2,tck=0.05,padj=1.4,cex.axis=0.9)
            }
        }
    }
    graphics::par(mar=c(0,0,0,0))
    graphics::plot(0,type="n",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
    graphics::text(1,0,paste(yvar,"~",xvar),cex=2)
}
