\name{ S3 plot method for transfer functions }
\alias{ plot.tfa }
\title{ Plot transfer function }
\description{
    Plot a transfer function
}
\usage{
    \method{plot}{tfa}(x, xvar=NULL, yvar=NULL, ...)
}
\arguments{
  \item{x}{
    an object of class 'tfa' created using \code{\link{tfa}} or \code{\link{inertia.tfa}}.
  }
  \item{xvar,yvar}{
    (optional) the variables to plot on the x and y axes.  May be \code{"p"}, \code{"lambda"}
    or \code{"inertia"}.  Defaults to \code{xvar="p"} and \code{yvar="lambda"} for objects created
    using \code{tfa} and \code{xvar="p"} and \code{yvar="inertia"} for objects created using
    \code{inertia.tfa}.
  }
  \item{...}{
    arguments to be passed to methods: see \code{\link{par}} and \code{\link{plot}}.
  }
}
\details{
    Plots transfer functions (class tfa) created using \code{\link{tfa}} or \code{\link{inertia.tfa}}.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    \code{\link{tfa}}, \code{\link{inertia.tfa}}
}
\examples{
    # Create a 3x3 matrix
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Calculate the transfer function of A[3,2] given a range of lambda
    evals<-eigen(A)$values
    lmax<-which.max(Re(evals))
    lambda<-Re(evals[lmax])
    lambdarange <- seq(lambda-0.1, lambda+0.1, 0.01)
    transfer<-tfa(A, d=c(0,0,1), e=c(0,1,0), lambdarange=lambdarange)
    transfer

    # Plot the transfer function (defaults to lambda~p)
    plot(transfer)

    # Create an initial stage structure
    initial <- c(1,3,2)
    initial

    # Calculate the transfer function of upper bound on inertia 
    # given a perturbation to A[3,2]
    transfer2<-inertia.tfa(A, d=c(0,0,1), e=c(0,1,0), bound="upper",
                          prange=seq(-0.6,0.4,0.01))
    transfer2

    # Plot the transfer function (defaults to p and inertia in
    # this case)
    plot(transfer2)

    # Plot lambda and inertia
    plot(transfer2, xvar="lambda", yvar="inertia")
}
\keyword{ plot }
\keyword{ transfer function }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection }
\keyword{ transient dynamics }
\keyword{ population projection matrix }
\keyword{ PPM }