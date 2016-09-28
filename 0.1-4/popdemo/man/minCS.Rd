\name{ minCS }
\alias{ minCS }
\title{ Calculate the mininum column sum of a matrix }
\description{
    Calculate the minimum column sum of a matrix.
}
\usage{
    minCS(A)
}
\arguments{
  \item{A}{
    a numeric matrix of any dimension.
  }
}
\details{
    \code{minCS} finds the smallest column sum of any numeric matrix.
}
\value{
    The minimum colmn sum of \code{A}.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\examples{
    # Create a 3x3 PPM
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Find the minimum column sum of A
    minCS(A)
}
\keyword{ minCS }
\keyword{ ecology }
\keyword{ demography }
\keyword{ transient dynamics }
\keyword{ population projection matrix }
\keyword{ PPM }