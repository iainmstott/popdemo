\name{ S3 plot method for matrices of transfer functions }
\alias{ plot.tfamatrix }
\title{ Plot matrices of transfer functions }
\description{
    Plot a matrix of transfer functions
}
\usage{
    \method{plot}{tfamatrix}(x, xvar=NULL, yvar=NULL, mar=c(1.1,1.1,0.1,0.1), ...)
}
\arguments{
  \item{x}{
    an object of class 'tfamatrix' created using \code{\link{tfamatrix}} or \code{\link{inertia.tfamatrix}}.
  }
  \item{xvar,yvar}{
    (optional) the variables to plot on the x and y axes.  May be \code{"p"}, \code{"lambda"}
    or \code{"inertia"}.  Defaults to \code{xvar="p"} and \code{yvar="lambda"} for objects created
    using \code{tfamatrix} and \code{xvar="p"} and \code{yvar="inertia"} for objects created using
    \code{inertia.tfamatrix}.
  }
  \item{mar}{
    the margin limits on the plots: see \code{\link{par}}
  }
  \item{...}{
    arguments to be passed to methods: see \code{\link{par}} and \code{\link{plot}}.
  }
}
\details{
    Plots matrices of transfer functions (class tfamatrix) created using 
    \code{\link{tfamatrix}} or \code{\link{inertia.tfamatrix}}.  The
    plot is laid out to correspond with the nonzero entries of the matrix
    used to generate the transfer functions, for easy visual comparison of
    how perturbation affects different matrix elements.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    \code{\link{tfamatrix}}, \code{\link{inertia.tfamatrix}}
}
\examples{
    # Create a 3x3 matrix
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Create an initial stage structure
    initial <- c(1,3,2)
    initial

    # Calculate the matrix of transfer functions using default arguments
    tfmat<-tfamatrix(A)
    tfmat

    # Plot the result
    plot(tfmat)

    # Calculate the matrix of transfer functions for inertia and 
    # specified initial stage structure using default arguments
    tfmat2<-inertia.tfamatrix(A,vector=initial)
    tfmat

    # Plot the result
    plot(tfmat2)

    # Plot inertia as a function of lambda
    plot(tfmat2, xvar="lambda", yvar="inertia")
}
\keyword{ plot }
\keyword{ transfer function }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection }
\keyword{ transient dynamics }
\keyword{ population projection matrix }
\keyword{ PPM }