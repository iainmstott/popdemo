################################################################################
#' Desert tortoise matrix
#'
#' @description
#' Population projection matrix for the desert tortoise \emph{Gopherus 
#' agassizii} with medium fecundity.
#'
#' @details 
#' The matrix is based on a population in the Western Mojave desert.
#' Stages are based on age and size (carapace length in mm):\cr
#' Yearling (age 0-1)\cr
#' Juvenile 1 (<60 mm)\cr
#' Juvenile 2 (90-99mm)\cr
#' Immature 1 (100-139mm)\cr
#' Immature 2 (140-179mm)\cr
#' Subadult (180-207mm)\cr
#' Adult 1 (208-239mm)\cr
#' Adult 2 (>240mm).
#'
#' @references
#' Doak et al. (1994) Ecol. Appl., 4, 446-460.
#'
#' @examples
#'   # read in data
#'   data(Tort)
#'   Tort
#'
#' @concept 
#' KEYWORDS
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
