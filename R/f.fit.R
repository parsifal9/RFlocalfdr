#' fit a spline to the histogram of imp
#'
#' 
#' @param imp the variable importances
#' @param df the degrees of freedom for the spline fit
#' @keywords spline
#' @export
#' @return a list with the following components
#' - "x" -- midpoints of the histogram
#' - "zh" -- a histogram object as returned by "hist"
#' - "f.spline" -- the spline fit. The fit is given by a glm mode glm(zh$counts ~ splines:::ns(x), poisson)
#' - "counts" the counts from the histogram
#' @md
#' @examples
#' data(imp1)                                                       
#' res <- f.fit(imp1)                                               
#' plot(res$zh, xlab="importances", main="histogram of importances")
#' points(res$midpoints,res$counts, col="grey90")                           
#' lines(res$zh$breaks[-1],res$f.spline,col="blue", lwd=3)          
#' legend(1.5,800,c("spline fit"), col="blue", lwd=3)               

f.fit <-function(imp,df=10,debug.flag=0,temp.dir=NULL){
    # do we want to set the number of breakpoints as a parameter
    # check diagnostic plot
    bre = 120
    lo <- min(imp)
    up <- max(imp)
    zzz <- pmax(pmin(imp, up), lo)
    breaks <- seq(lo, up, length = bre)
    zh <- hist(zzz, breaks = breaks, plot = FALSE)
    x <- (breaks[-1] + breaks[-length(breaks)])/2 #midpoints
    y <- zh$counts  
    f.spline <- glm(y ~ splines:::ns(x, df = df), poisson)$fit

    if (debug.flag > 0){
        png(paste(temp.dir,"/histogram_of_variable_importances.png",sep=""))
        hist.of.data<-hist(zzz, breaks = breaks,freq=TRUE, xlab="importance",axes=FALSE,
                           main="histogram of variable importances" )
        lines(x,f.spline,col="red",lwd="4")
        lines(x,hist.of.data$counts,col="green",lwd="4")
        legend("topright",c("spline","counts"),col=c("red","green"),lty=1,lwd=4)
        box()
        dev.off()
    }
    
    temp <- list(x, y, zh, f.spline)
    names(temp) <- c("midpoints", "counts", "zh", "f.spline")
    temp
}
