################################################################################
#' Plot transfer function
#'
#' @description
#' Plot a transfer function
#'
#' @param x an object of class 'tfa' (transfer function analysis) created using 
#' \code{\link{tfa_lambda}} or \code{\link{tfa_inertia}}.
#' @param xvar,yvar (optional) the variables to plot on the x and y axes. May
#' be \code{"p"}, \code{"lambda"} or \code{"inertia"}. Defaults to
#' \code{xvar="p"} and \code{yvar="lambda"} for objects created using
#' \code{tfa_lambda} and \code{xvar="p"} and \code{yvar="inertia"} for 
#' objects created using \code{tfa_inertia}.
#' @param ...  arguments to be passed to methods: see \code{\link{par}} and
#' \code{\link{plot}}.
#'
#' @details 
#' \code{plot.tfa} plots transfer functions (class \code{tfa}) created using 
#' \code{\link{tfa_lambda}} or \code{\link{tfa_inertia}}.
#' 
#'
#' @seealso
#' Constructor functions: \code{\link{tfa_lambda}}, \code{\link{tfa_inertia}}
#'
#' @examples
#'   # Create a 3x3 matrix
#'   ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#'
#'   # Calculate the transfer function of A[3,2] given a range of lambda
#'   evals <- eigen(A)$values
#'   lmax <- which.max(Re(evals))
#'   lambda <- Re(evals[lmax])
#'   lambdarange <- seq(lambda-0.1, lambda+0.1, 0.01)
#'   ( transfer <- tfa_lambda(A, d=c(0,0,1), e=c(0,1,0), lambdarange=lambdarange) )
#'
#'   # Plot the transfer function
#'   plot(transfer)
#'
#'   # Create an initial stage structure
#'   ( initial <- c(1,3,2) )
#' 
#'   # Calculate the transfer function of upper bound on inertia 
#'   # given a perturbation to A[3,2]
#'   ( transfer<-tfa_inertia(A, d=c(0,0,1), e=c(0,1,0), bound="upper",
#'                           prange=seq(-0.6,0.4,0.01)) )
#' 
#'   # Plot the transfer function (defaults to inertia ~ p)
#'   plot(transfer)
#'
#'   # Plot inertia against lambda
#'   plot(transfer, xvar="lambda", yvar="inertia")
#' 
#' @concept transfer function
#' @concept systems control
#' @concept nonlinear
#' @concept perturbation
#' @concept population viability
#' @concept PVA
#' @concept ecology
#' @concept demography
#' 
#' @method plot tfa
#' @export
plot.tfa <- function(x,xvar=NULL,yvar=NULL,...){
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
        yv<-x[[length(x)]]
        yvar<-names(x)[length(x)]
    }
    else{
        if(yvar=="p") yv<-x$p
        if(yvar=="lambda") yv<-x$lambda
        if(yvar=="inertia") yv<-x$inertia
    }
    graphics::plot(xv,yv,type="l",xlab=xvar,ylab=yvar,...)
}
