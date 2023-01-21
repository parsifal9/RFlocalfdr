#' local fdr 
#'
#' calculate the local 
#' @param f  f 
#' @param x  x 
#' @param FUN FUN
#' @keywords cats
#' @export
#' @examples
#' local.fdr()
local.fdr <-function(f,x,FUN=dburr, p0= 1, debug.flag=0,plot.string="",temp.dir=NULL, ...){
    f0<-FUN(x, ...)
    f<-f$f.spline
    f <- (sum(f0) * f)/sum(f)
    fdr <- pmin((p0 * f0)/f,1)

    if (debug.flag >1 ){
        png(paste(temp.dir,"/local_fdr_",plot.string,".png",sep=""))
        plot(fdr)
        abline(h=0.2)
        dev.off()
    }
    
    fdr
}


