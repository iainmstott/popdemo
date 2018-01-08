################################################################################
#' Polar bear matrices
#' 
#' Matrix projection model for the polar bear \emph{Ursus maritimus}, with 
#' 5 matrices corresponding to years 2001-2005. The matrices are based on a 
#' population in the southern Beaufort Sea. Stages are based on age and 
#' reproductive status:\cr
#' Stage-1: 2-year-old\cr
#' Stage 2: 3-year-old\cr
#' Stage 3: 4-year-old\cr
#' Stage 4: adult (5+ years old), available to breed\cr
#' Stage 5: adult, with cub (0-1 years old)\cr
#' Stage 6: adult, with yearling (1-2 years old).\cr\cr
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
#' Hunter et al. (2010) Ecology, 91, 2883-2897.
#' 
#' @examples
#'   #read in data
#'   data(Pbear)
#'   Pbear
#' 
#' @concepts 
#' data polar bear matrix Ursus maritimus
#'
Pbear <- list(matrix(c(0,0,0,0,0,0.5811,
                       0.9858,0,0,0,0,0,
                       0,0.9858,0,0,0,0,
                       0,0,0.9858,0.5061,0.3791,0.9918,
                       0,0,0,0.4857,0.0681,0,
                       0,0,0,0,0.5433,0),
                     ncol=6,byrow=TRUE),
              matrix(c(0,0,0,0,0,0.58,
                       0.9842,0,0,0,0,0,
                       0,0.9842,0,0,0,0,
                       0,0,0.9842,0.4563,0.3654,0.9911,
                       0,0,0,0.5348,0.0808,0,
                       0,0,0,0,0.5435,0),
                     ncol=6,byrow=TRUE),
              matrix(c(0,0,0,0,0,0.5379,
                       0.9415,0,0,0,0,0,
                       0,0.9415,0,0,0,0,
                       0,0,0.9415,0.284,0.3081,0.9662,
                       0,0,0,0.6822,0.1384,0,
                       0,0,0,0,0.516,0),
                     ncol=6,byrow=TRUE),
              matrix(c(0,0,0,0,0,0.2773,
                       0.6578,0,0,0,0,0,
                       0,0.6578,0,0,0,0,
                       0,0,0.6578,0.5367,0.4243,0.7587,
                       0,0,0,0.222,0.0327,0,
                       0,0,0,0,0.2689,0),
                     ncol=6,byrow=TRUE),
              matrix(c(0,0,0,0,0,0.3165,
                       0.7034,0,0,0,0,0,
                       0,0.7034,0,0,0,0,
                       0,0,0.7034,0.7225,0.4328,0.7943,
                       0,0,0,0.0718,0.0081,0,
                       0,0,0,0,0.3254,0),
                     ncol=6,byrow=TRUE))
names(Pbear) <- c("Y2001","Y2002","Y2003","Y2004","Y2005")
dimnames(Pbear[[1]]) <- list(c("2YO", "3YO", "4YO", "Ad", "AdC", "AdY"),
                             c("2YO", "3YO", "4YO", "Ad", "AdC", "AdY"))
dimnames(Pbear[[2]]) <- list(c("2YO", "3YO", "4YO", "Ad", "AdC", "AdY"),
                             c("2YO", "3YO", "4YO", "Ad", "AdC", "AdY"))
dimnames(Pbear[[3]]) <- list(c("2YO", "3YO", "4YO", "Ad", "AdC", "AdY"),
                             c("2YO", "3YO", "4YO", "Ad", "AdC", "AdY"))
dimnames(Pbear[[4]]) <- list(c("2YO", "3YO", "4YO", "Ad", "AdC", "AdY"),
                             c("2YO", "3YO", "4YO", "Ad", "AdC", "AdY"))
dimnames(Pbear[[5]]) <- list(c("2YO", "3YO", "4YO", "Ad", "AdC", "AdY"),
                             c("2YO", "3YO", "4YO", "Ad", "AdC", "AdY"))
