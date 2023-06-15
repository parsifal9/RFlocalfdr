#' my_ranger_PIMP
#' based on the PIMP function from the vita package. ‘PIMP’ implements the test approach of Altmann et al. (2010) for
#' the permutation variable importance measure ‘VarImp’ returned by the randomForest package (Liaw and Wiener (2002)) for
#' classification and regression.
#'
#' my_PIMP applies the same method as PIMP but to the MDI (mean decrease in impurity)  variable importance
#' (mean decrease in Gini index for classification and  mean decrease in MSE for regression).
#' my_ranger_PIMP applies the same method to the ranger RF package
#' @param  X
#' @param  y
#' @param  rForest
#' @param  S
#' @param  parallel
#' @param  ncores
#' @param  seed
#' @md
#' @export
#' @return an object of class PIMP
#' @examples
#' library(RFlocalfdr)
#' library(vita) #vita: Variable Importance Testing Approaches
#' data(smoking)
#' ?smoking 
#' y<-smoking$y
#' y<-factor(y)
#' smoking_data<-smoking$rma
#'
#' cl.ranger <- ranger(y=y, x=smoking_data,mtry = 3,num.trees = 1000, importance = 'impurity')
#' system.time(pimp.varImp.cl<-my_ranger_PIMP(X,y,cl.ranger,S=10, parallel=TRUE, ncores=3)) 
#' pimp.t.cl <- PimpTest(pimp.varImp.cl,para = FALSE)
#' aa <- summary(pimp.t.cl,pless = 0.05)
#' length(which(aa$cmat2[,"p-value"]< 0.05))
#' hist(aa$cmat2[,"p-value"],breaks=20)

my_ranger_PIMP <- function (X, y, rForest, S = 100, parallel = FALSE, ncores = 0, 
    seed = 123, ...) 
{
    if (!inherits(rForest, "ranger")) 
        stop("rForest is not of class ranger")
    mtry <- rForest$mtry
    num.trees <- rForest$num.trees
    n <- nrow(X)
    p <- ncol(X)
    if (parallel) {
        d_ncores = parallel::detectCores()
        if (ncores > d_ncores) 
            stop("ncores: The number of cores is too large")
        if (ncores == 0) {
            ncores = max(c(1, floor(d_ncores/2)))
        }
        if ("L'Ecuyer-CMRG" != RNGkind()[1]) {
            cat("\n The random number generator was set to L'Ecuyer-CMRG !! \n")
            RNGkind("L'Ecuyer-CMRG")
        }
    }
    else {
        ncores = 1
    }
    set.seed(seed)
    classRF = is.factor(y)
    if (classRF) {
        y.s = replicate(S, sample(as.integer(y)))
        i = 1:S
        varip <- parallel::mclapply(i, function(i) {
           aa<- ranger::ranger(y=as.factor(y.s[, i]), x=X,
                mtry = mtry, importance = "impurity", num.trees = num.trees
                )$variable.importance
#           aa<- ranger::ranger(as.factor(y.s[, i])~., X,
#                mtry = mtry, importance = "impurity", num.trees = num.trees
#                )$variable.importance
        }, mc.cores = ncores)
        varip <- simplify2array(varip)
    }
     else {
         y.s = replicate(S, sample(y))
         i = 1:S
         varip = parallel::mclapply(i, function(i) {
             ranger::ranger(y=y.s[, i],x= X, mtry = mtry, 
                 importance = "impurity", num.trees = num.trees)$variable.importance
#             ranger::ranger(y.s[, i]~., X, mtry = mtry, 
#                 importance = "impurity", num.trees = num.trees)$variable.importance
         }, mc.cores = ncores)
         varip = simplify2array(varip)
     }
    ##dimNames <-  names(rForest$variable.importance)

    ## if (classRF) {
    ##     VarImp <- rForest$variable.importance
    ##     }
    ##     else {
    ##         VarImp = rForest$importance[, 1]
    ##         }
            
    out = list(VarImp = matrix(rForest$variable.importance, ncol = 1, 
                               dimnames = list(names(rForest$variable.importance))), PerVarImp = varip, 
               type = if (classRF) {
                          "classification"
                      } else {
                          "regression"
                      }, call = match.call())
    class(out) = "PIMP"
    return(out)
}