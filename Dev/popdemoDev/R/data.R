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
#' \itemize{
#'  \item Doak et al. (1994) Ecol. Appl., 4, 446-460.
#' }
#'
#' @examples
#'   # read in data and view
#'   data(Tort)
#'   Tort
#'
#' @concept 
#' data desert tortoise matrix Gopherus agassizzii
#'
"Tort"

################################################################################
#' Polar bear matrices
#' 
#' Matrix projection model for the polar bear \emph{Ursus maritimus}, with 
#' 5 matrices corresponding to years 2001-2005. The matrices are based on a 
#' population in the southern Beaufort Sea. During 2001-2003, ice conditions were
#' classified as "good", but in 2004-2005, ice conditions were classified as 
#' "poor". Poor ice conditions lead to worse population performance. Stages are 
#' based on age and 
#' reproductive status:\cr
#' Stage-1: 2-year-old\cr
#' Stage 2: 3-year-old\cr
#' Stage 3: 4-year-old\cr
#' Stage 4: adult (5+ years old), available to breed\cr
#' Stage 5: adult, with cub (0-1 years old)\cr
#' Stage 6: adult, with yearling (1-2 years old).
#' 
#' The population structure is 
#' \code{c(0.106, 0.068, 0.106, 0.461, 0.151, 0.108))}
#' 
#' @docType data
#' 
#' @usage data(Pbear)
#' 
#' @format
#' List object containing matrices.
#' 
#' @references
#' \itemize{
#'  \item Hunter et al. (2010) Ecology, 91, 2883-2897.
#' }
#' 
#' @examples
#'   #read in data
#'   data(Pbear)
#'   Pbear
#' 
#' @concept
#' data polar bear matrix Ursus maritimus
#'
"Pbear"
