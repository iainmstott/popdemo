#' @name Stochastic-Matrix-Classes
#' @aliases CompadreDDM CompadreMSM
#' 
#' @title Stochastic Matrix Classes in popdemo
#' 
#' @description Stochastic matrix classes in \code{Rcompadre} and \code{popdemo} are broken
#' down into two types: stochastic by matrix selection and stochastic by element
#' selection. Within the latter category, there are three different types supported:
#' density dependent, environmentally dependent, and their combination (note that
#' only density dependence is implemented currently, but the last two are coming
#' in future versions). 
#' 
#' Stochastic by matrix selection matrices are of class \code{CompadreMSM} and
#' only contain a list of matrices filled with actual numbers. Additional information
#' for each one should be supplied in the call to \code{project}.
#' 
#' Stochastic by element selection matrices are comprised of two lists - one for
#' constants and the population vector (the \code{dataList}), and one for 
#' expressions that calculate vital rates based on the population vector, as well
#' as the expression for the matrix itself (the \code{matExprs} and 
#' \code{matrixExpr}, respectively). The population vector is optional and can be
#' set to \code{NULL}. Options for stage-biased vectors and random draws from a
#' dirichlet distribution are provided in \code{project}.
#'  
#' Entries in \code{matExprs} that are not the \code{matrixExpression} should read
#' like R code. If a value in one expression is a function of other constants, 
#' \code{makeMatExprs} will recognize this and generate the correct order for
#' evaluating the expressions. For example, fecundity of the second stage 
#' class could density dependent and defined as \deqn{f_2 = e^{bf_0 + bf_1 * u}}
#' This should be translated into an expression as \code{f_2 = exp(bf_0 + bf_1 * u)}
#' (A poisson model with slope \code{bf_1} and intercept \code{bf_0}). \code{u} is
#' the sum of all non-dormant stages in the population vector. Thus, if you supply
#' and additional expression, \code{u = sum(stage1, stage2, stage3)}, then the 
#' expression for \code{u} will get evaluated first and substituted into the 
#' expression for \code{f_2}. More complicated expressions (provided they are valid
#' R code) can also be supplied as long as they return a single numeric value. 
#' You do not need to specify these in any particular 
#' order in calls to \code{makeMatExprs}. This is most useful after constructing
#' a model for a vital rate - generating an expression that calls 
#' \code{predict(model,...)} is valid.
#' 
#' See the vignette on writing density 
#' dependent functions for additional examples (\code{browseVignettes('popdemo')}).
#'  
#' The matrix expression can be filled with a
#' mixture of symbols (e.g. \code{s_2 * f}) and numeric values, provided that
#' the symbols are defined in either \code{matExprs} or named in the 
#' \code{dataList}. It should be supplied in the form of \code{c(a11, a12, a13, ...)} 
#' and in row major order. Entries can be numbers or expressions provided the 
#' expression is valid R code. \code{\link{makeMatExprs}} is provided to assist with 
#' getting the format correct.
#' 
#' @examples
#' # Create a CompadreMSM object from a list of matrices 
#' 
#' data(Pbear)
#' my_msm <- CompadreMSM(A = Pbear)
#' 
#' # If there is only one matrix, it will need to be put into a list
#' 
#' data(Tort)
#' tort_msm <- CompadreMSM(A = list(Tort))
#' 
#' # Creating a DDM requires a set of matrix expressions and a data list.
#' # This example makes use of Pardini et al (2009) Complex dynamics and control of
#' # invasive bienniel Alliaria petiolata (garlic mustard). Ecologogical Applications
#' 
#' my_matExprs <- makeMatExprs(
#'   s_2 = 1/(1 + exp(bs2_2 * u_i + bs2_1 * t_i + bs2_0)),
#'   s_3 = exp(bs3_1 * log(r + 1)),
#'   f = exp(bf_1 * a + bf_0),
#'   u_i = r * a, # density dependent terms. These appear in s_2, s_3, and f
#'   t_i = r + a, # r and a are the second and third stages in initPopVector
#'   matrixExpr =
#'       c(
#'         1 - g_2, 0, v * (1-g_1) * f,
#'         g_2 * s_1, 0, v * g_1 * s_1 * f,
#'         0, s_2 * s_3, 0
#'       ),
#'    matrixDimension = 3
#' )
#' 
#' my_dataList <- makeDataList(
#'   v = 0.8228,
#'   g_1 = 0.5503,
#'   g_2 = 0.3171,
#'   bs2_2 = 0.0016,
#'   bs2_1 = -0.0664,
#'   bs2_0 = -0.156,
#'   bs3_1 = -0.289,
#'   bf_1 = -0.0389,
#'   bf_0 = 7.489,
#'   s_1 = 0.5,
#'   initPopVector = c(s = 10, r = 0, a = 0)
#' )
#' 
#' my_DDM <- CompadreDDM(dataList = my_dataList,
#'                       matExprs = my_matExprs)
#' 
#' 
NULL

#' @rdname Stochastic-Matrix-Classes
#' @export CompadreDDM
#' @slot dataList list of named constant values for substitution in matrix. For 
#' additional examples, see \code{\link{makeDataList}}. 
#' @slot matExprs list of expressions that define the density dependent components.
#' This can be constructed by hand, but it is safer to construct this using 
#' \code{makeMatExprs}.
#' For additional examples, see \code{\link{makeMatExprs}}.
#' 

CompadreDDM <- setClass('CompadreDDM',
                        slots = c(dataList = 'list',
                                  matExprs = 'list'))

ValidCompadreDDM <- function(object) {
  

  if(!'matrixExpr' %in% names(object@matExprs)) {
    stop('User-defined "CompadreDDM" objects must contain a "matrixExpr" expression.\n',
         'See ?makeMatExprs for more details.')
  }

  # need to add ability to check for potential negative matrix values.
  TRUE
}

setValidity("CompadreDDM", ValidCompadreDDM)

#' @rdname Stochastic-Matrix-Classes
#' @export CompadreMSM
#' @slot A A list of matrices

CompadreMSM <- setClass('CompadreMSM', 
                        slots = c(A = 'list'))

validCompadreMSM <- function(object) {
  anyNegative <- function(x) any(x < 0)
  mats <- object@A
  
  negs <- vapply(mats, anyNegative, logical(1))
  squares <- vapply(mats, .isSquare, logical(1))
  # Can add more here
  
  errors <- character()
  if(any(negs)) {
    negInd <- which(negs)
    errors <- c(errors, paste('Matrices ', paste(negInd, collapse = ', '), 
                              ' contain negative entries.', sep = ""))
  }
  
  if(any(!squares)) {
    sqInd <- which(!squares)
    errors <- c(errors, paste('Matrices ', paste(sqInd, collapse = ', '), 
                              ' are not square.', sep = ""))
  }
  
  ifelse(length(errors) > 0,
         stop(.errorConstructor(errors), call. = FALSE),
         TRUE)
}



setValidity('CompadreMSM', validCompadreMSM)

# CompadreESM is not yet ready, but this code is! uncomment this section once
# we're ready to write a project(CompadreESM) method
# #' @rdname Stochastic-Matrix-Classes
# #' 
# #' @inheritParams CompadreDDM
# #' @slot envExprs additional expressions specifying environmental effects
# #' @export 

# CompadreESM <- setClass('CompadreESM',
#                         slots = c(data_list = 'list',
#                                   matExprs = 'list',
#                                   envExprs = 'list'))
# 
# 
# ValidCompadreESM <- function(object) {
# 
#   if(!'matrixExpr' %in% names(object@matExprs)) {
#     stop('User-defined "CompadreESM" objects must contain a "matrixExpr" expression.\n',
#          'See ?makeMatExprs for more details.')
#   }
#   
#   # need to add ability to check for potential negative matrix values.
#   TRUE
# }
# 
# setValidity("CompadreESM", ValidCompadreESM)

#' @rdname Stochastic-Matrix-Classes
#' @param object An object of class \code{CompadreMSM} or \code{CompadreDDM}.
#' @export

setGeneric("matrixExpr", function(object) {
  standardGeneric('matrixExpr')
})

#' @rdname Stochastic-Matrix-Classes
#' @export
setMethod('matrixExpr', signature = c("CompadreDDM"),
          function(object) {
            temp <- rlang::quo_get_expr(object@matExprs$matrixExpr)
            # out <- expr_to_char_mat(temp) # need to write this
            return(temp)
          }
)

# #' @rdname Stochastic-Matrix-Classes
# #' @export
# setMethod('matrixExpr', signature = c("CompadreESM"),
#           function(object) {
#             temp <- rlang::quo_get_expr(object@matExprs$matrixExpr)
#             # out <- expr_to_char_mat(temp) # need to write this
#             return(temp)
#           }
# )

#' @rdname Stochastic-Matrix-Classes
#' @export
setGeneric('matFuns', function(object) {
  standardGeneric('matFuns')
})

#' @rdname Stochastic-Matrix-Classes
#' @export
setMethod('matFuns', signature = c("CompadreDDM"),
          function(object) {
            ind <- which(names(object@matExprs) == 'matrixExpr')
            funs <- object@matExprs[-ind]
            out <- lapply(funs, rlang::quo_get_expr)
            return(out)
          }
)


# #' @rdname Stochastic-Matrix-Classes
# #' @export
# setMethod('matFuns', signature = c("CompadreESM"),
#           function(object) {
#             ind <- which(names(object@matExprs) == 'matrixExpr' |
#                            names(object@matExprs) == 'matDim')
#             funs <- object@matExprs[-ind]
#             out <- lapply(funs, rlang::quo_get_expr)
#             return(out)
#           }
# )


