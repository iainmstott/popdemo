\name{ S3 plot method for projections }
\alias{ plot.projection }
\title{ Plot population dynamics }
\description{
    Plot dynamics of a specified population 
    projection matrix (PPM) model.
}
\usage{
    \method{plot}{projection}(x, labs=TRUE, ...)
}
\arguments{
  \item{x}{
    an object of class 'projection' created using \code{\link{project}}.
  }
  \item{labs}{
    logical: if \code{FALSE}, then lines are not labelled in stage-biased projections.
  }
  \item{...}{
    arguments to be passed to methods: see \code{\link{par}} and \code{\link{plot}}.
  }
}
\details{
    Plots population dynamics (time series of density) for objects of class 'projection'
    created using \code{\link{project}}.  The method is particularly useful for sets of stage-biased
    projections (which \code{project} creates by default when \code{vector="n"}).\cr
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    \code{\link{project}}
}
\examples{
    # Create a 3x3 PPM
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Create an initial stage structure
    initial <- c(1,3,2)
    initial

    # plot stage-biased dynamics of A over 70 intervals using
    # standardised dynamics and log y axis
    plot(project(A,time=70,standard.A=TRUE,standard.vec=TRUE),
         log="y",ylab="Years",xlab="Density")

    # plot a projection of a specified initial stage structure
    # for 10 intervals
    plot(project(A,vector=initial,time=50))

}
\keyword{ project }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection }
\keyword{ transient dynamics }
\keyword{ population projection matrix }
\keyword{ PPM }