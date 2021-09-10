#' local fdr 
#'
#' This function allows you to express your love of cats.
#' @param f  f 
#' @param FUN FUN
#' @keywords cats
#' @export
#' @examples
#' local.fdr()
local.fdr <-function(f,FUN=dburr, ...){
    f0<-FUN(x, ...)
    f<-f$f.spline
    
    p0<-1
    f <- (sum(f0) * f)/sum(f)
    fdr <- pmin((p0 * f0)/f,1)
    plot(fdr)
    abline(h=0.2)
    fdr
}


