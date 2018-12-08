
#' @rdname Density-Dependent-Constructor-Funs
#' @aliases makeMatExprs makeDataList
#' @title Helpers to construct CompadreDDMs
#' 
#' @description These functions assist with creation of \code{CompadreDDM}s by ensuring that
#' each expression and associated values are evaluated in the correct order.
#' 
#' 
#' @param ... For \code{makeDataList}, named constant values that are 
#' substituted into the \code{matExprs} expressions. These can be regression 
#' coefficients, values for fixed vital rates, or any other value that appears
#' in \code{matExprs}. For \code{makeMatExprs}, this is a set of named expressions.
#' The left hand side of each one should be a parameter that appears either 
#' in the \code{matrixExpr} or in another expression for a density dependent vital
#' rate.
#' @param initPopVector optionally, an initial population vector to start a
#' \code{Projection} with. Values in this should be named with values corresponding
#' to the \code{matrixExpr} or \code{matExprs}. See details for warnings on not 
#' specifying this!
#' 
#' @return \code{makeDataList} returns a named list with parameter values.
#' \code{makeMatExpr} returns a named list of expressions that are used
#' to calculate values of matrix elements at each iteration, and the matrix itself.
#' For those who are curious or wish to construct these by hand, these expressions 
#' are stored as \code{\link[rlang]{quos}} with \code{env} slot set to \code{empty}.
#' The environment is reassigned to a specific evaluation environment during
#' iteration and the user-specified one (if it is specified at all) will be 
#' overridden. 
#' 
#' @examples 
#' # This example makes use of Pardini et al (2009) Complex dynamics and control of
#' # invasive bienniel Alliaria petiolata (garlic mustard). Ecologogical Applications
#' 
#' # makeMatExprs can take raw expressions, even when the parameters in them are
#' # defined by other expressions. This allows you to focus on specifying each
#' # density dependent model correctly without worrying about the specifics of
#' # function evaluation.
#' 
#' exprs <- makeMatExprs(
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
#' # Use this for constants and the initial population vector. All values
#' # here should appear in the expressions above (either in a vital rate or 
#' # matrix element). Notice that the initial popopulation vector is named such
#' # that stages listed there match the values in "u_i" and "t_i".
#' 
#' data <- makeDataList(
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
#' alliariaDDM <- CompadreDDM(dataList = data,
#'                            matExprs = exprs)
#' 
#' project(alliariaDDM, time = 100)
#' 
#' # Now, one without a population vector. In this case, we subtract the 
#' # initPopVector entry from our dataList and 
#' 
#' vec_less_data <- data[-(length(data))]
#' vec_less_exprs <- makeMatExprs(
#'   s_2 = 1/(1 + exp(bs2_2 * u_i + bs2_1 * t_i + bs2_0)),
#'   s_3 = exp(bs3_1 * log(V2 + 1)), # notice that this is now switched from "r" to V2
#'   f = exp(bf_1 * V3 + bf_0),      # and this is now V3 instead of a
#'   u_i = V2 * V3,                  # V2 * V3 instead of r * a
#'   t_i = V2 + V3,                  # V2 + V3 instead of r + a
#'   matrixExpr =
#'     c(
#'       1 - g_2, 0, v * (1-g_1) * f,
#'       g_2 * s_1, 0, v * g_1 * s_1 * f,
#'       0, s_2 * s_3, 0
#'     ),
#'   matrixDimension = 3
#' )
#' 
#' vecless_alliariaDDM <- CompadreDDM(dataList = vec_less_data,
#'                                    matExprs = vec_less_exprs)
#' \dontrun{
#' # Test out dirichlet vectors and make sure plotting works
#' vecless_dirichlet_alliaria <- project(vecless_alliariaDDM,
#'                                       vector = 'diri', 
#'                                       draws = 500, 
#'                                       alpha.draws = 'unif')
#' }
#' 
#' @export

makeDataList <- function(...,
                         initPopVector = NULL) {
  
  data <- list(...)
  data$popVec <- initPopVector
  return(data)
}


#' @rdname Density-Dependent-Constructor-Funs
#' 
#' @inheritParams makeDataList
#' 
#' @param matrixExpr An expression that specifies the form of the density dependent
#' matrix. Supplied as c(...) and in row-major format (e.g. forms correctly when
#' \code{byrow = TRUE}).
#' @param matrixDimension An integer that specifies how many rows and columns 
#' are in the matrix.
#' 
#' @details Note that if you do not plan to pass an initial population vector
#' to \code{makeDataList}, then the names of the vector will be automatically 
#' generated as \code{V1}, \code{V2}, etc. This means that the expressions in
#' \code{makeMatExprs} that are functions of the population vectors \strong{must
#' use these names}.
#' 
#' @importFrom rlang enquos enquo quo_is_null
#' @export

makeMatExprs <- function(...,
                         matrixExpr = NULL,
                         matrixDimension = NULL) {
  
  matExprs <- rlang::enquos(...)
  matExprs$matrixExpr <- rlang::enquo(matrixExpr)
  
  if(rlang::quo_is_null(matExprs$matrixExpr)) {
    stop('Must supply an expression for "matrixExpr".',
         call. = FALSE)
  }
  
  if(is.null(matrixDimension)) {
    stop('Must supply an integer for "matrixDimension".')
  }
  
  matExprs$matrixExpr <- .wrapMatrixCall(matExprs$matrixExpr, 
                                         matrixDimension)
  
  textList <- .wrapEvalTidys(matExprs)
  
  out <- lapply(textList, .wrapQuos)
  
  out$matDim <- matrixDimension
  return(out)
}

#' @noRd
#' @importFrom rlang quo_text enquo parse_expr
.wrapMatrixCall <- function(matrixExpr, matrixDimension) {
  
  # get text, wrap it in the matrix call with correct row count
  matexpr <- rlang::quo_text(matrixExpr)
  
  text_call <- paste('matrix(', 
                     matexpr,
                     ', nrow = ', 
                     matrixDimension,
                     ', byrow = TRUE)', sep = "")
  
  # back to expression!
  parsed_call <- rlang::parse_expr(text_call)
  
  # enquo it so we can save it for later!
  new_quo <- rlang::enquo(parsed_call)
  
  return(new_quo)
}

#' @noRd
#' @importFrom rlang quo_text eval_tidy
.wrapEvalTidys <- function(matExprs) {
  
  # turn quos into strings, extract the LHS variable names
  textList <- lapply(matExprs, rlang::quo_text)
  LHS <- names(textList)
  
  for(i in seq_along(textList)) {
    for(j in seq_along(LHS)){
      
      # substitute eval_tidy(var) if anything in right hand side appears in
      # left hand side
      regExpression <- .makeMatRegexpr(LHS[j])
      
      textList[[i]] <- gsub(regExpression,
                             paste0('rlang::eval_tidy(', LHS[j], ')'),
                             textList[[i]])
      
    }
  }
  return(textList)
}

#' @noRd
.makeMatRegexpr <- function(var) {
  paste0('(\\b', var, '\\b)')
}

#' @noRd
#' @importFrom rlang quo
.wrapQuos <- function(matExprsWEvals) {
  
  string <- paste0('rlang::quo(', matExprsWEvals, ')')
  out <- rlang::parse_expr(string)
  enquo(out)
  
}

#' @noRd
#' @importFrom rlang env env_bind env_bind_lazy !!! quo_set_env

.initMatEnv <- function(matExprs, dataList) {
  # insulated environment for the iterations to take place
  evalEnv <- rlang::env()
  
  # remove the matrixDimension slot, it's not a quosure! 
  matDimInd <- which(names(matExprs) == 'matDim')
  matExprs <- matExprs[-matDimInd]
  
  # Set this as the environment for each quosure in matExprs
  matExprs <- lapply(matExprs, function(x) quo_set_env(x, evalEnv))
  
  # add data to evalEnv
  rlang::env_bind(evalEnv,
                  !!! dataList)
  
  # create lazy bindings to the evaluation environment
  rlang::env_bind_lazy(evalEnv,
                       !!! matExprs,
                       .eval_env = evalEnv)
  return(evalEnv)
}

.checkCurrentMat <- function(currentMat) {
  
  if(any(currentMat < 0)) {
    stop('Some element of the matrix in the current iteration',
         'is negative.\n',
         'Check expressions supplied to matExprs in the',
         'CompadreDDM object.')
  }
  
  # Would be good to add check for colSums > 1 on a U matrix, but not sure how
  # to split out the fecundity values.
  
}

#' @noRd
.getStageNames <- function(popVec) {
  # can only be a matrix or a vector if it is user supplied
  # This is so wildly inelegant - there must be a better way...
  if(is.matrix(popVec)) { 
    if(!is.null(dimnames(popVec)[[1]])) {
      out <- dimnames(popVec)[[1]]
    } else {
      out <- paste('V', 1:dim(x)[[1]])
    }
  } else {
    
    if(!is.null(names(popVec))) {
      out <- names(popVec)
    } else {
      out <- paste('V', 1:length(popVec))
    }
  }
  
  return(out)
}

.initDDPopVec <- function(vector, 
                          time,
                          matDim,
                          stageNames,
                          draws,
                          alpha.draws,
                          standard.vec) {
  
  if(is.numeric(vector)) {
    vectorSwitch <- 'user'
  } else {
    vectorSwitch <- vector
  }
  
  popVec <- switch(vectorSwitch,
                   'diri' = .initDiriVector(matDim = matDim,
                                            time = time,
                                            stageNames = stageNames,
                                            draws = draws,
                                            alpha.draws = alpha.draws),
                   'n' = .initBiasVector(matDim = matDim,
                                         time = time,
                                         stageNames = stageNames),
                   'user' = .initUserVector(vector = vector,
                                            time = time,
                                            matDim = matDim, 
                                            stageNames = stageNames,
                                            standard.vec = standard.vec),
                   # not yet implemented, waiting to see what
                   # the database version of these looks like
                   # 'db' = .initDbVector(vector, 
                   #                      matDim, 
                   #                      stageNames)
                   )
  
  if(is.numeric(vector)) {
    if((is.array(popVec) & dim(popVec)[3] == 1)) vecType <- 'single'
    if(dim(popVec)[3] > 1) vecType <- 'multiple'
  } else {
    vecType <- vector
  }
  
  projType <- 'stochastic'
  
  out <- list(popVec = popVec,
              projType = projType,
              vecType = vecType)
  return(out)
  
}
