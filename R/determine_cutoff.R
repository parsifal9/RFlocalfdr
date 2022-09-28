#' count the number of times each variable is used in the a ranger forest
#'
#' count the number of times each varaible is used in the a ranger forest.
#' help(treeInfo) warns
#' "splitvarID -- ID of the splitting variable, 0-indexed. Caution, the variable order changes if the formula interface is used"
#' this should be investigated
#' @param  imp importances
#' @param  t2 number of times each variable is used in the a ranger forest.
#' @param  cutoff
#' @param  Q 
#' @param  plot
#' @export
#' @examples
#' library(ranger)
#' rf1 <-ranger(Species ~ ., data = iris,importance="impurity",seed=123)
#' imp<-log(importance(rf1))
#' t2 <- count_variables(rf1)
#' aa <- determine_cutoff(imp, t2, cutoff=c(0,1,2), Q = 0.75, plot = NULL)
determine_cutoff <- function(imp, t2, cutoff=c(0,1,4,10,15,20), Q = 0.75, plot = NULL){
    res1  <- matrix(0, length(cutoff), 3)
    steps <- cutoff
    old.par <- par()
    par(mfrow=c(2,2))    # set the plotting area into a 2*2 array

    for ( ii in 1:length(steps )  ) {
        cat("i=", ii, "cutoff =", steps[ii], "\n")
        temp <- imp[t2 > steps[ii]]
        temp <- temp[temp != -Inf]
        temp <- temp - min(temp) + .Machine$double.eps
        f_fit <- f.fit(temp )
        y <- f_fit$zh$density
        x <- f_fit$x

        C <- quantile(temp,probs=Q)
        df2 <- data.frame(x[x < C], y[x < C])

        initial.estimates <- fit.to.data.set.wrapper(df2, temp)
        if(class(initial.estimates) != "nls"){
            return(res1)
        }
        initial.estimates <- summary(initial.estimates)$parameters
        initial.estimates <- data.frame( initial.estimates)
        t1<-my.dsn(x, xi = initial.estimates$Estimate[ 1], 
                   omega = initial.estimates$Estimate[ 2],  lambda =  initial.estimates$Estimate[ 3]
                   )
        t1<-t1[x < C]

        if (!is.null(plot)){
            if (steps[ii] %in% plot){
                plot(density(temp),main=paste("cutoff = ",steps[ii]))
                lines(df2$x, df2$y,col="red",lwd=2)
                abline(v=C,col="red" )
                lines(df2$x,t1,col="blue",lwd=2)
            }
        }
        res1[ii,1] <- sum((df2$y-t1)^2) 
        res1[ii,2] <- sum(abs(df2$y-t1))
        res1[ii,3] <- max(abs(df2$y-t1)) 
    }
     par(old.par)
    res1
}
