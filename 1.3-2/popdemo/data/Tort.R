################################################################################
#' Desert tortoise matrix
#'
#' Matrix Projection Model for the desert tortoise \emph{Gopherus agassizii} 
#' with medium fecundity. The matrix is based on a population in the Western 
#' Mojave desert. Stages are based on age and size 
#' (carapace length in mm):\cr
#' Stage 1: Yearling (age 0-1)\cr
#' Stage 2: Juvenile 1 (<60 mm)\cr
#' Stage 3: Juvenile 2 (90-99mm)\cr
#' Stage 4: Immature 1 (100-139mm)\cr
#' Stage 5: Immature 2 (140-179mm)\cr
#' Stage 6: Subadult (180-207mm)\cr
#' Stage 7: Adult 1 (208-239mm)\cr
#' Stage 8: Adult 2 (>240mm).
#'
#' @docType data
#'
#' @usage data(Tort)
#'
#' @format
#' Object of class \code{matrix}
#' 
#' @references
#' Doak et al. (1994) Ecol. Appl., 4, 446-460.
#'
#' @examples
#'   # read in data and view
#'   data(Tort)
#'   Tort
#'
#' @concept 
#' data desert tortoise matrix Gopherus agassizzii
#'
Tort <- matrix(c(0,0,0,0,0,1.3,1.98,2.57,
                 0.716,0.567,0,0,0,0,0,0,
                 0,0.149,0.567,0,0,0,0,0,
                 0,0,0.149,0.604,0,0,0,0,
                 0,0,0,0.235,0.56,0,0,0,
                 0,0,0,0,0.225,0.678,0,0,
                 0,0,0,0,0,0.249,0.851,0,
                 0,0,0,0,0,0,0.016,0.86),
                 ncol=8,byrow=TRUE)
dimnames(Tort) <- list(c("Yr", "J1", "J2", "I1", "I2", "SA", "A1", "A2"),
                       c("Yr", "J1", "J2", "I1", "I2", "SA", "A1", "A2"))
