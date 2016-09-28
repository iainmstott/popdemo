\name{ Tort }
\alias{ Tort }
\title{ Desert tortoise matrix }
\description{
    Population projection matrix for the desert tortoise Gopherus agassizii
    with medium fecundity.
}
\usage{
    data(Tort)
}
\details{
    The matrix is based on a population in the Western Mojave desert.
    Stages are based on age and size (carapace length in mm):\cr
    Yearling (age 0-1)\cr
    Juvenile 1 (<60 mm)\cr
    Juvenile 2 (90-99mm)\cr
    Immature 1 (100-139mm)\cr
    Immature 2 (140-179mm)\cr
    Subadult (180-207mm)\cr
    Adult 1 (208-239mm)\cr
    Adult 2 (>240mm).
}
\references{
    Doak et al. (1994) Ecol. Appl., 4, 446-460.
}
\examples{
    # read in data
    data(Tort)
    Tort
}
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection matrix }
\keyword{ PPM }