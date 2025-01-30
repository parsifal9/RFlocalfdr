#' evaluate a measure that can be used to determining a significance level for the Mean Decrease in Impurity measure returned by a Random Forest model
#'
#' @param  imp vector of MDI variable importances 
#' @param  t2 number of times each variable is used in the a ranger forest. Returned by             
#'         count_variables for a ranger RF
#' @param  cutoff values to evaluate 
#' @param  Q -- we examine the fit up to the quartile Q
#' @param  plot for 4 selected values of the cutoff. The plot contains
#' - The data (black) and  the  fitted density  (red)
#' - The skew-normal fit (blue)
#' - The quantile Q (vertical red line)
#' @param verbose verbose=0, no output to screen
#'                verbose=1, track the cutoff value being used
#' @param try.counter passed to fit.to.data.set.wrapper
#' @md
#' @return res a matrix if size length(cutoff) by 3.
#' We model the histogram of imp with a kernel density estimate, y.
#' Let t1 be  fitted values of the skew normal. Then res contains three columns
#' - sum((y-t1)^2)
#' - sum(abs(y-t1)) and
#' - max(abs(y-t1)),
#' evaluated up to the quantile Q
#' @export
#' @examples
#' \dontrun{
#' data(imp20000)
#' imp <- log(imp20000$importances)
#' t2<- imp20000$counts
#' length(imp)
#' hist(imp,col=6,lwd=2,breaks=100,main="histogram of importances")
#' res.temp <- determine_cutoff(imp, t2 ,cutoff=c(0,1,2,3),plot=c(0,1,2,3),Q=0.75,try.counter=1)
#' 
#' plot(c(0,1,2,3),res.temp[,3])
#' # the minimum is at 1 so
#' imp<-imp[t2 > 1]
#' 
#' qq <- plotQ(imp,debug.flag = 0)
#' ppp<-run.it.importances(qq,imp,debug=0)
#' aa<-significant.genes(ppp,imp,cutoff=0.2,debug.flag=0,do.plot=2, use_95_q=TRUE)
#' length(aa$probabilities) #11#
#' names(aa$probabilities)
#' #[1] "X101"   "X102"   "X103"   "X104"   "X105"   "X2994"  "X9365"  "X10718"
#' # [9] "X13371" "X15517" "X16460"
#' # so the observed FDR is 0.54
#' }
#' 
#' \dontrun{
#' library(ranger)
#' library(RFlocalfdr.data)
#' data(smoking)
#' y<-smoking$y
#' smoking_data<-smoking$rma
#' y.numeric <-ifelse((y=="never-smoked"),0,1)
#' rf1 <- ranger::ranger(y=y.numeric ,x=smoking_data,importance="impurity",seed=123, num.trees = 10000,     
#'          classification=TRUE)
#' t2 <-count_variables(rf1)
#' imp<-log(rf1$variable.importance)
#' plot(density((imp)))
#' # Detemine a cutoff to get a unimodal density.
#' res.temp <- determine_cutoff(imp, t2 ,cutoff=c(1,2,3,4),plot=c(1,2,3,4),Q=0.75)
#' plot(c(1,2,3,4),res.temp[,3])
#' }

determine_cutoff <- function (imp, t2, cutoff = c(0, 1, 4, 10, 15, 20), Q = 0.75, 
    plot = NULL, verbose = 0, try.counter=3)
{
    res1 <- matrix(0, length(cutoff), 3)
    steps <- cutoff
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(2, 2))
    for (ii in 1:length(steps)) {
        if (verbose > 0) {
            message("i = ", ii, ", cutoff = ", steps[ii], "\n")
        }
        temp <- imp[t2 > steps[ii]]
        temp <- temp[temp != -Inf]
        temp <- temp - min(temp) + .Machine$double.eps
        f_fit <- f.fit(temp)
        y <- f_fit$zh$density
        x <- f_fit$midpoints
        C <- quantile(temp, probs = Q)
        df2 <- data.frame(x[x < C], y[x < C])
        initial.estimates <- fit.to.data.set.wrapper(df2, temp, try.counter = try.counter)
        if (!inherits(initial.estimates, "try-error")){
            initial.estimates <- summary(initial.estimates)$parameters
            initial.estimates <- data.frame(initial.estimates)
            t1 <- my.dsn(x, xi = initial.estimates$Estimate[1], omega = initial.estimates$Estimate[2], 
                         lambda = initial.estimates$Estimate[3])
            t1 <- t1[x < C]
        } else {
            t1 <- NULL
           }
        if (!is.null(plot) && inherits(t1, "numeric")) {
            if (steps[ii] %in% plot) {
                plot(density(temp), main = paste("cutoff = ", steps[ii]))
                lines(df2$x, df2$y, col = "red", lwd = 2)
                abline(v = C, col = "red")
                lines(df2$x, t1, col = "blue", lwd = 2)
            }
        }
        if (inherits(t1, "numeric")){
            res1[ii, 1] <- sum((df2$y - t1)^2)
            res1[ii, 2] <- sum(abs(df2$y - t1))
            res1[ii, 3] <- max(abs(df2$y - t1))
        } else {
            res1[ii, 1] <- NaN
            res1[ii, 2] <- NaN
            res1[ii, 3] <- NaN
           }

    }
    res1
}

