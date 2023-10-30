#' count the number of times each variable is used in a ranger random forest
#'
#' * count the number of times each variable is used in a ranger random forest.
#' * help(treeInfo) warns "splitvarID -- ID of the splitting variable, 0-indexed.
#'    Caution, the variable order changes if the formula interface is used"
#' However this should be investigated
#' @param object a ranger forest object
#' @returns a table (0-indexed) giving the number of times each variable was used in the random forest 
#' @keywords counts
#' @export
#' @examples
#' library(ranger)
#' rf1 <- ranger(Species ~ ., data = iris,importance="impurity", seed=123)
#' count_variables(rf1)
#' rf2 <- ranger(dependent.variable.name = "Species", data = iris,seed=123)
#' count_variables(rf2)
#' rf3<- ranger(y = iris[, 5], x = iris[, -5],seed=123)
#' count_variables(rf3)
#' 

count_variables<-function(object) {
  if (!inherits(object, "ranger")) {
    stop("Error: Invalid class of input object.")
  }
  aa1<- unlist(object$forest$split.varIDs)
  aa2<-unlist(lapply(object$forest$child.nodeIDs,f <- function(x) { x[[1]] }))
  t1<- table(aa1[aa2 != 0]) 
  temp<-object$num.independent.variables
  t2<- table(factor(aa1[aa2 != 0],levels=0:(temp-1)))
  t2
}

