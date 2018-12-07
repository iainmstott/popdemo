#' @name Stochastic-Matrix-Classes
#' @aliases CompadreDDM CompadreMSM
#' 
#' @title Stochastic Matrix Classes 
#' 
#' @description Stochastic matrix classes in \code{Rcompadre} and \code{popdemo} can be broken
#' down into two types: stochastic by matrix selection and stochastic by element
#' selection. Within the latter category, there are three different types supported:
#' density dependent, environmentally dependent, and their combination (note that
#' only density dependence is implemented currently, but the last two are coming
#' in future versions). 
#' 
#' Stochastic by matrix selection matrices are of class \code{CompadreMSM} and
#' only contain a list of matrices filled with actual numbers. Additional information
#' for each one should be supplied in call to \code{project}.
#' 
#' Stochastic by element selection matrices are comprised of two lists - one for
#' constants and the population vector (the \code{dataList}), and one for 
#' expressions that calculate vital rates based on the population vector, as well
#' as the expression for the matrix itself (the \code{matExprs} and 
#' \code{matrixExpr}, respectively). The population vector is optional and can be
#' set to \code{NA}. Options for stage-biased vectors and random draws from a
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
#' expression for \code{s_2}. More complicated expressions (provided they are valid
#' R code) can also be supplied as long as they return a scalar numeric value. 
#' You do not need to specify these in any particular 
#' order in calls to \code{makeMatExprs}. See the vignette on writing \code{matExprs}
#' functions for additional examples (\code{browseVignettes('popdemo')}).
#'  
#' The matrix expression can be filled with a
#' mixture of symbols (e.g. \code{s_2 * f}) and numeric values, provided that
#' the symbols are defined in either \code{matExprs} or named in the 
#' \code{dataList}. See the examples in \code{\link{makeMatExprs}}.
#' 
NULL

#' @rdname Stochastic-Matrix-Classes
#' @export
#' @slot data_list list of named constant values for substitution in matrix. For 
#' additional examples, see \code{\link{makeDataList}}. 
#' @slot matExprs list of expressions that define the density dependent components.
#' This can be constructed by hand, but it is safer to construct this using 
#' \code{makeMatExprs}.
#' For additional examples, see \code{\link{makeMatExprs}}.
#' 

CompadreDDM <- setClass('CompadreDDM', contains = "list",
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
#' @export
#' @slot A a list of matrices

CompadreMSM <- setClass('CompadreMSM', contains = 'list',
                        slots = c(A = 'list'))

#' @rdname Stochastic-Matrix-Classes
#' 
#' @inheritParams CompadreDDM
#' @slot env_exprs additional expressions specifying environmental effects
#' @export 

CompadreESM <- setClass('CompadreESM', contains = "list",
                        slots = c(data_list = 'list',
                                  matExprs = 'list',
                                  envExprs = 'list'))


ValidCompadreESM <- function(object) {

  if(!'matrixExpr' %in% names(object@matExprs)) {
    stop('User-defined "CompadreESM" objects must contain a "matrixExpr" expression.\n',
         'See ?makeMatExprs for more details.')
  }
  
  # need to add ability to check for potential negative matrix values.
  TRUE
}

setValidity("CompadreESM", ValidCompadreESM)

#' @rdname Stochastic-Matrix-Classes
#' 
#' 

#' @rdname Stochastic-Matrix-Classes
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

#' @rdname Stochastic-Matrix-Classes
#' @export
setMethod('matrixExpr', signature = c("CompadreESM"),
          function(object) {
            temp <- rlang::quo_get_expr(object@matExprs$matrixExpr)
            # out <- expr_to_char_mat(temp) # need to write this
            return(temp)
          }
)

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


#' @rdname Stochastic-Matrix-Classes
#' @export
setMethod('matFuns', signature = c("CompadreESM"),
          function(object) {
            ind <- which(names(object@matExprs) == 'matrixExpr')
            funs <- object@matExprs[-ind]
            out <- lapply(funs, rlang::quo_get_expr)
            return(out)
          }
)


