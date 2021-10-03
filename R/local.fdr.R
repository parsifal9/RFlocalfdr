#' local fdr 
#'
#' This function allows you to express your love of cats.
#' @param f  f 
#' @param x  x 
#' @param FUN FUN
#' @keywords cats
#' @export
#' @examples
#' local.fdr()
local.fdr <-function(f,x,FUN=dburr, p0= 1, ...){
    f0<-FUN(x, ...)
    f<-f$f.spline
    f <- (sum(f0) * f)/sum(f)
    fdr <- pmin((p0 * f0)/f,1)
    plot(fdr)
    abline(h=0.2)
    fdr
}


