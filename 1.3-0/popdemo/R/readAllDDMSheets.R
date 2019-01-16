#' Reads all the sheets in Compadre Density-dependent matrix Excel sheet
#' 
#' @description this function reads all the pages in an Excel sheet
#' and returns a list of data frames. The length of the list is equal to the
#' number of the sheets in the file. Unfortunately, sheet names are not
#' preserved so you'll need to insert those manually if you want to use them.
#' @param file the file path to the Excel sheet.
#' 
#' @return A list of data frames
#' 
#' @importFrom readxl read_excel
#' @export

readAllSheets <- function(file) {
  
  sheets <- readxl::excel_sheets(file)
  
  out <- lapply(sheets,
                FUN = function(x) readxl::read_excel(file, sheet = x))
  
  names(out) <- sheets
  
  # padEnv <- rlang::new_environment(data = out)
  
  return(out)
  
}


#' Converts a list of data frames into a CompadreDDM
#' 
#' @param parameters A list produced by \code{readAllSheets}.
#' 
#' @return an object of class \code{CompadreDDM}.
#' 
#' @export 

list2DDM <- function(parameters) {

  CompadreDDM(dataList = .makeTabDataList(parameters[[1]]),
              matExprs = .makeTabExprList(parameters[[2]]))
  
}
 
# extracts constants and creates a named list for a ComapdreDDM
.makeTabDataList <- function(dataList) {
  
  out <- list()
  for(i in seq_len(dim(dataList)[1])) {
    out[[i]] <- eval(rlang::parse_expr(dataList$value[i]))
    names(out)[i] <- dataList$parameter[i]
  }
  
  return(out)

}

# Extracts the expressions and parses them into the right format
.makeTabExprList <- function(matExprs) {
  
  dotsInd <- which(!matExprs$expr_name %in% c('matrixDimension', 'matrixExpr'))
  matExprInd <- which(matExprs$expr_name %in% 'matrixExpr')
  matDimInd <- which(matExprs$expr_name %in% 'matrixDimension')
  
  # Get dots for makeMatExprs
  dots <- rlang::parse_exprs(matExprs$expr[dotsInd])
  names(dots) <- matExprs$expr_name[dotsInd]
  
  # get matrixExpr
  matExpr <- rlang::parse_expr(matExprs$expr[matExprInd])

  # get integer for matrixDimension
  matDim <- eval(rlang::parse_expr(matExprs$expr[matDimInd]))
  
  out <- makeMatExprs(!!! dots, matrixExpr = !! matExpr, matrixDimension = matDim)
  
  return(out)
}

