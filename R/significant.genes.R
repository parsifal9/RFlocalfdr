#' significant.genes
#'
#' This function allows you to express your love of cats.
#' @param object object 
#' @param imp importances
#' @param cutoff cutoff
#' @param do.plot do.plot 
#' @param debug.flag debug.flag
#' @keywords cats
#' @export
#' @examples
#' cat_function()
significant.genes <- function(object,imp,cutoff=0.2,do.plot=0,debug.flag=0){
    
    start.x<- which.min(abs(object$fdr-mean(imp)))
    ww<-which.min(abs(object$fdr[start.x:119]-cutoff)) 
    num.sig.genes <-sum(imp> object$x[as.numeric(names(ww))])

     if (do.plot == 1){
         plot(object$x,object$fdr)
         abline(h=cutoff)
         abline(v=object$x[as.numeric(names(ww))])
     }
    
    if(debug.flag==1){
        cat(sum(imp> object$x[as.numeric(names(ww))]),"sum(imp1> x[as.numeric(names(ww))])","\n")
    }
     if (do.plot==2){
         hist(imp, breaks = 200,freq=FALSE)
         abline(v=sn::qsn(0.95,  xi=object$estimates[1], omega=object$estimates[2], alpha= object$estimates[3]))
         abline(v=object$x[as.numeric(names(ww))])
     }
     if(debug.flag==2){
         cat(names(imp1)[imp1> object$x[as.numeric(names(ww))]],"\n")
     }
    
    a1<- match(names(imp1)[imp1> object$x[as.numeric(names(ww))]],names(imp1))
    ppp<- 1-sn::psn(imp1[a1], xi=object$estimates[1], omega=object$estimates[2], alpha= object$estimates[3])
    names(ppp) <-names(imp1)[imp1> object$x[as.numeric(names(ww))]]
     if (debug.flag==2){
         cat(length(ppp),"\n") 
     }
    ppp
}
