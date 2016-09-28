\name{ R2Matlab }
\alias{ R2Matlab }
\title{ Convert matrices into Matlab style strings }
\description{
    Convert \R objects of class matrix into character
    strings that represent the matrix in a Matlab style
}
\usage{
    R2Matlab(A, noquote=FALSE)
}
\arguments{
  \item{A}{
    a numeric matrix of any dimension
  }
  \item{noquote}{
    (optional) if \code{noquote=TRUE} then the returned character vector is printed without quotes.
  }
}
\details{
    Matlab reads matrices using a unique one-line notation that can prove useful for storage in databases
    and importing multiple matrices into a program at once, amongst other applications.  This notation 
    is by row, with "[" and "]" to specify the beginning and end of the matrix respectively, ";" to 
    specify a new row and a space between each matrix element.  Thus, the \R matrix created using 
    \code{matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)} may is equivalent to
    [0 1 2;0.5 0.1 0;0 0.6 0.6].\cr

    \code{R2Matlab} takes an \R object of class matrix and uses string-building techniques to convert
    it into a Matlab-style character string that may be useful for exporting into databases.
}
\value{
    Object of class character representing \code{A} in a Matlab style.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    \code{\link{Matlab2R}}
}
\examples{
    # Create a 3x3 PPM
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Code the matrix in a Matlab style
    R2Matlab(A)

    # Print without quotes
    R2Matlab(A, noquote=TRUE)
}
\keyword{ R2Matlab }
\keyword{ Matlab }
\keyword{ matrix }