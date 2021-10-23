#' fit a spline 
#'
#' fit a spline 
#' @param zz the varaible importances
#' @param df the degrees of freedom for the spline fit
#' @keywords spline
#' @export
#' @examples

f.fit <-function(zz,df=10,debug.flag=0,temp.dir=NULL){
    bre = 120
    lo <- min(zz)
    up <- max(zz)
    zzz <- pmax(pmin(zz, up), lo)
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
    
    hist.of.data<-hist(zzz, breaks = breaks, plot = FALSE)
    f.hist<-hist.of.data$counts
    temp <- list(x, zh, f.spline,hist.of.data$counts)
    names(temp) <- c("x", "zh", "f.spline", "counts")
    temp
}
