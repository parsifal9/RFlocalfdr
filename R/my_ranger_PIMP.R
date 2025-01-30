#' my_ranger_PIMP
#' based on the PIMP function from the vita package. ‘PIMP’ implements the test approach of Altmann et al. (2010) for
#' the permutation variable importance measure ‘VarImp’ returned by the randomForest package (Liaw and Wiener (2002)) for
#' classification and regression.
#'
#' my_PIMP applies the same method as PIMP but to the MDI (mean decrease in impurity)  variable importance
#' (mean decrease in Gini index for classification and  mean decrease in MSE for regression).
#' my_ranger_PIMP applies the same method to the ranger RF package
#' @param  X, data matrix of size n by p
#' @param  y, class labels for classification (factor) or real values for regression. Of length n
#' @param  rForest, an object of class ranger, importance must be set to "impurity".
#' @param  S, The number of permutations for the response vector ‘y’. Default is ‘S=100
#' @param  parallel Should the PIMP-algorithm run parallel?  Default is
#'          ‘parallel=FALSE’ and the number of cores is set to one. The
#'          parallelized version of the PIMP-algorithm are based on
#'            mclapply and so is not available on Windows
#' @param  ncores, The number of cores to use, i.e. at most how many child
#'           processes will be run simultaneously. Must be at least one,
#'           and parallelization requires at least two cores.  If
#'           ‘ncores=0’, then the half of CPU cores on the current host
#'           are used.
#' @param  seed a single integer value to specify seeds. The "combined
#'           multiple-recursive generator" from L'Ecuyer (1999) is set as
#'           random number generator for the parallelized version of the
#'           PIMP-algorithm.  Default is ‘ seed = 123’.
#' @param ... additional arguments passed to ranger
#' @md
#' @export
#' @return an object of class PIMP
#' @examples
#' \dontrun{ 
#' library(RFlocalfdr.data)
#' library(ranger)
#' library(vita) #vita: Variable Importance Testing Approaches
#' data(smoking)
#' ?smoking 
#' y<-smoking$y
#' y<-factor(y)
#' smoking_data<-smoking$rma
#'
#' cl.ranger <- ranger::ranger(y=y, x=smoking_data,mtry = 3,num.trees = 1000, importance = 'impurity')
#' system.time(pimp.varImp.cl<-my_ranger_PIMP(smoking_data,y,cl.ranger,S=10, parallel=TRUE, ncores=2))
#' #CRAN limits the number of cores available to packages to 2, for performance reasons.
#' pimp.t.cl <- vita::PimpTest(pimp.varImp.cl,para = FALSE)
#' aa <- summary(pimp.t.cl,pless = 0.05)
#' length(which(aa$cmat2[,"p-value"]< 0.05))
#' hist(aa$cmat2[,"p-value"],breaks=20)
#' }

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
            message("\n The random number generator was set to L'Ecuyer-CMRG !! \n")
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
