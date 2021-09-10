#' fit a spline 
#'
#' fit a spline 
#' @param zz the varaible importances
#' @param df the degrees of freedom for the spline fit
#' @keywords spline
#' @export
#' @examples
f.spline <-function(zz,df=10){
    bre = 120
    lo <- min(zz)
    up <- max(zz)
    zzz <- pmax(pmin(zz, up), lo)
    breaks <- seq(lo, up, length = bre)
    zh <- hist(zzz, breaks = breaks, plot = F)
    x <- (breaks[-1] + breaks[-length(breaks)])/2 #midpoints
    y <- zh$counts  
    f <- glm(y ~ splines:::ns(x, df = df), poisson)$fit
    hist.of.data<-hist(zzz, breaks = breaks,freq=TRUE, xlab="importance",axes=FALSE,
                       main="histogram of VI" )
    lines(x,f,col="red",lwd="4")
    box()
    
    f.spline<-f
    
    f.hist<-hist.of.data$counts
    lines(x,f.spline,col="blue",lwd="4")
    lines(x,f.hist,col="green",lwd="4")
    temp<-list(x,zh,f.spline)
    names(temp)<-c("x","zh","f.spline")
    temp
}
