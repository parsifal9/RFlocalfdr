#' my_PIMP
#' based on the PIMP function from the vita package. ‘PIMP’ implements the test approach of Altmann et al. (2010) for
#' the permutation variable importance measure ‘VarImp’ in a random
#'     forest for classification and regression.
#'
#' @param  X
#' @param  y
#' @param  rForest
#' @param  S
#' @param  parallel
#' @param  ncores
#' @param  seed
#' @md
#' @export
#' @return res a matrix if size length(cutoff) by 3.
#' We model the histogram of imp with a kernel density estimate, y.
#' Let t1 be  fitted values of the skew normal. Then res contains three columns
#' - sum((y-t1)^2)
#' - sum(abs(y-t1)) and
#' - max(abs(y-t1)),
#' evaluated up to the quantile Q
#' @examples
#' data(ch22)                                                                                 
#' ? ch22                                                                                     
#' #document how the data set is created                                                      
#' plot(density(log(ch22$imp)))                                                               
#' t2 <-ch22$C                                                                                
#' imp<-log(ch22$imp)                                                                         
#' #Detemine a cutoff to get a unimodal density.                                              
#' res.temp <- determine_cutoff(imp, t2 ,cutoff=c(1,10,20,30),plot=c(1,10,20,30),Q=0.75)      
#' plot(c(1,2,3,4),res.temp[,3])                                                              
#'                                                                                            
#' res.temp <- determine_cutoff(imp, t2 ,cutoff=c(25,30,35,40),plot=c(25,30,35,40),Q=0.75)    
#' plot(c(25,30,35,40),res.temp[,3])                                                          
#' imp<-imp[t2 > 30]
#' qq <- plotQ(imp,debug.flag = 0)
#' ppp<-run.it.importances(qq,imp,debug=0)                                                       
#' aa<-significant.genes(ppp,imp,cutoff=0.2,debug.flag=0,do.plot=2)                           
#' length(aa$probabilities) # 6650                                                            
#' aa<-significant.genes(ppp,imp,cutoff=0.05,debug.flag=0,do.plot=2)                          
#' length(aa$probabilities) # 3653    

my_PIMP <- function(X, y, rForest, S = 100, parallel = FALSE, ncores = 0,
    seed = 123, ...)
{
    if (!inherits(rForest, "randomForest")) 
        stop("rForest is not of class randomForest")
    mtry <- rForest$mtry
    ntree <- rForest$ntree
    n <- nrow(X)
    p <- ncol(X)
    if (Sys.info()[["sysname"]] == "Windows" & parallel) {
        cat("\n The parallelized version of the PIMP-algorithm are not available on Windows !! \n")
        parallel = FALSE
    }
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
        if (rForest$type == "regression") 
            stop("rForest$type = regression !! y a factor ")
        y.s = replicate(S, sample(as.integer(y)))
        i = 1:S
        varip = parallel::mclapply(i, function(i) {
            randomForest::randomForest(X, as.factor(y.s[, i]), 
                mtry = mtry, importance = TRUE, ntree = ntree, 
                ...)[[9]][, 3]
        }, mc.cores = ncores)
        varip = simplify2array(varip)
    }
    else {
        if (rForest$type == "classification") 
            stop("rForest$type = classification !! y not a factor ")
        y.s = replicate(S, sample(y))
        i = 1:S
        varip = parallel::mclapply(i, function(i) {
            randomForest::randomForest(X, y.s[, i], mtry = mtry, 
                importance = TRUE, ntree = ntree, ...)[[7]][, 
                1]
        }, mc.cores = ncores)
        varip = simplify2array(varip)
    }
    dimNames = dimnames(rForest$importance)[[1]]

    if (classRF) {
        VarImp <- rForest$importance[, 3]
        }
        else {
            VarImp = rForest$importance[, 1]
            }
            
            out = list(VarImp = matrix(VarImp, ncol = 1, 
        dimnames = list(dimNames, "VarImp")), PerVarImp = varip, 
        type = if (classRF) {
            "classification"
        } else {
            "regression"
        }, call = match.call())
    class(out) = "PIMP"
    return(out)
}





my.ranger.PIMP <- function (X, y, rForest, S = 100, parallel = FALSE, ncores = 0, 
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
