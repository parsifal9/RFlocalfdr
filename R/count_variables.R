#' count the number of times each varaible is used in the a ranger forest
#'
#' count the number of times each varaible is used in the a ranger forest.
#' help(treeInfo) warns
#' "splitvarID -- ID of the splitting variable, 0-indexed. Caution, the variable order changes if the formula interface is used"
#' this should be investigated
#' @param object a ranger forest object
#' @keywords counts
#' @export
#' @examples
#' library(ranger)
#' rf1<-ranger(Species ~ ., data = iris,importance="impurity",seed=123)
#' count_variables(rf1)
count_variables<-function(object){
    if (!inherits(object, "ranger")) {
        stop("Error: Invalid class of input object.")
    }
    aa1<- unlist(object$forest$split.varIDs)
    aa2<-unlist(lapply(object$forest$child.nodeIDs,f<-function(x){x[[1]]}))
    t1<- table(aa1[aa2 != 0]) 
    temp<-object$num.independent.variables
    t2<- table(factor(aa1[aa2 != 0],levels=0:(temp-1)))
    t2
}






## library(ranger)
## rf1<-ranger(Species ~ ., data = iris,importance="impurity",seed=123)
## rf2 <- ranger(dependent.variable.name = "Species", data = iris,importance="impurity",seed=123)
## rf3<-  ranger(y = iris[, 5], x = iris[, -5],data = iris,importance="impurity",seed=123)

## aa1<- unlist(rf1$forest$split.varIDs)
## table(aa1)

## aa <- treeInfo(rf1, tree = 1)
## aa$splitvarID
## aa$splitvarName
## table(aa$splitvarID,aa$splitvarName)
## # so 0 1 2 3
## #Sepal.Length Sepal.Width Petal.Length  Petal.Width
## rf1$variable.importance
## #Sepal.Length  Sepal.Width Petal.Length  Petal.Width 
## #   10.101955     2.456635    43.068485    43.57366

## # help(treeInfo) splitvarID -- ID of the splitting variable, 0-indexed. Caution, the variable order changes if the formula interface is used. 


## aa1<- unlist(rf2$forest$split.varIDs)
## table(aa1)
## aa <- treeInfo(rf2, tree = 1)
## aa$splitvarID
## aa$splitvarName
## table(aa$splitvarID,aa$splitvarName)


## aa1<- unlist(rf3$forest$split.varIDs)
## table(aa1)
## rf1$forest$splitvarName
## aa <- treeInfo(rf3, tree = 1)
## aa$splitvarID
## aa$splitvarName
## table(aa$splitvarID,aa$splitvarName)



## object<-rf1
## aa1<- unlist(object$forest$split.varIDs)
## aa2<-unlist(lapply(object$forest$child.nodeIDs,f<-function(x){x[[1]]}))
## t1<- table(aa1[aa2 != 0]) 
## temp<-object$num.independent.variables
## t2<- table(factor(aa1[aa2 != 0],levels=1:temp))
## t2


## aa <- treeInfo(rf1, tree = 1)
## aa$splitvarID
## aa$splitvarName


## count_variables(rf1)
