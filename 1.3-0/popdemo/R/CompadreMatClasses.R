#' @name Stochastic-Matrix-Classes
#' 
#' @title Stochastic Matrix Classes 
#' 
#' @description Stochastic matrix classes in \code{Rcompadre} and \code{popdemo} can be broken
#' down into two type: stochastic by matrix selection and stochastic by element
#' selection. Within the latter category, there are three different types supported:
#' density dependent, environmentally dependent, and their combination. 
#' 
#' Stochastic by matrix selection matrices are of class \code{CompadreMSM} and
#' must contain a list of matrices. They can also contain additional information
#' such as initial population vectors, pre-defined sequences of integers that 
#' define the order in which projection matrices are selected, time intervals
#' for projection 
#' 
NULL

#' @rdname Stochastic-Matrix-Classes
#' @export
#' @slot data_list list of named constant values for substitution in matrix
#' @slot mat_exprs list of expressions that define the density dependent components
#' 

CompadreDDM <- setClass('CompadreDDM', contains = "list",
                        slots = c(dataList = 'list',
                                  matExprs = 'list'))

ValidCompadreDDM <- function(object) {
  

  if(!'matExpr' %in% names(object@matExprs)) {
    stop('User-defined "CompadreDDM" objects must contain a "matExpr" expression.\n',
         'See ?makeMatExprs for more details.')
  }
  
  # need to add ability to check for potential negative matrix values.
  TRUE
}

setValidity("CompadreDDM", ValidCompadreDDM)

#' @rdname Stochastic-Matrix-Classes
#' @export
#' @slot A a list of matrices
#' @slot vector optional initial population vector or matrix of vectors
#' @slot time optional integer specifying default number of iterations for
#' \code{project}
#' @slot standard.A optional logical specifying whether to standardize
#' matrices in \code{A} by their respective dominant eigenvalues
#' @slot standard.vec optional logical specifying whether to  standardize population
#' vector to sum to 1
#' @slot Aseq optional integer sequence specifying order to choose matrices in
#' \code{project}
#' @slot Astart optional integer or character specifying position or name of 
#' matrix in A to use as starting point for projections.
#' @slot draws number of draws to use if \code{vector = 'diri'}
#' @slot alpha.draws blah blah blah
#' @slot PREcheck blah blah blah
#' 

CompadreMSM <- setClass('CompadreMSM', contains = 'list',
                        slots = c(A = 'list',
                                  vector = NULL,
                                  time = NULL,
                                  standard.A = NULL,
                                  standard.vec = NULL,
                                  Aseq = NULL,
                                  Astart = NULL,
                                  draws = NULL,
                                  alpha.draws = NULL,
                                  PREcheck = NULL))

#' @rdname Stochastic-Matrix-Classes
#' 
#' @inheritParams CompadreDDM
#' @slot env_exprs additional expressions specifying environmental effects
#' @export 

CompadreEPM <- setClass('CompadreEPM', contains = "list",
                        slots = c(data_list = 'list',
                                  mat_exprs = 'list',
                                  env_exprs = 'list'))


ValidCompadreEPM <- function(object) {

  if(!'mat_expr' %in% names(object@mat_exprs)) {
    stop('User-defined "CompadreEPM" objects must contain a "mat_expr" expression.\n',
         'See ?make_mat_exprs for more details.')
  }
  
  # need to add ability to check for potential negative matrix values.
  TRUE
}

setValidity("CompadreEPM", ValidCompadreEPM)

#' @rdname Stochastic-Matrix-Classes
#' 
#' 

#' @rdname Stochastic-Matrix-Classes
#' @export
setGeneric("mat_expr", function(object) {
  standardGeneric('mat_expr')
})

#' @rdname Stochastic-Matrix-Classes
#' @export
setMethod('mat_expr', signature = c("CompadreDDM"),
          function(object) {
            temp <- rlang::quo_get_expr(object@mat_exprs$mat_expr)
            # out <- expr_to_char_mat(temp) # need to write this
            return(temp)
          }
)

#' @rdname Stochastic-Matrix-Classes
#' @export
setMethod('mat_expr', signature = c("CompadreEPM"),
          function(object) {
            temp <- rlang::quo_get_expr(object@mat_exprs$mat_expr)
            # out <- expr_to_char_mat(temp) # need to write this
            return(temp)
          }
)

#' @rdname Stochastic-Matrix-Classes
#' @export
setGeneric('mat_funs', function(object) {
  standardGeneric('mat_funs')
})

#' @rdname Stochastic-Matrix-Classes
#' @export
setMethod('mat_funs', signature = c("CompadreDDM"),
          function(object) {
            ind <- which(names(object@mat_exprs) == 'mat_expr')
            funs <- object@mat_exprs[-ind]
            out <- lapply(funs, rlang::quo_get_expr)
            return(out)
          }
)


#' @rdname Stochastic-Matrix-Classes
#' @export
setMethod('mat_funs', signature = c("CompadreEPM"),
          function(object) {
            ind <- which(names(object@mat_exprs) == 'mat_expr')
            funs <- object@mat_exprs[-ind]
            out <- lapply(funs, rlang::quo_get_expr)
            return(out)
          }
)


