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

significant.genes <-function(object,imp,cutoff=0.2,do.plot=0,debug.flag=0){
    start.x <- which.min(abs(object$x-mean(imp)))
    ww<-which.min(abs(object$fdr[start.x:119]-cutoff)) 
    num.sig.genes <-sum(imp> object$x[as.numeric(names(ww))])

    if (do.plot == 1){
        hist(imp,col=6,lwd=2,breaks=100,main="",freq=FALSE,xlab="importances",ylab="density",ylim=c(0,1))
         lines(object$x,object$fdr)
         abline(h=cutoff)
         abline(v=object$x[as.numeric(names(ww))])
     }
    
    if(debug.flag==1){
        cat(sum(imp> object$x[as.numeric(names(ww))]),"sum(imp> x[as.numeric(names(ww))])","\n")
    }
     if (do.plot==2){
         hist(imp,col=6,lwd=2,breaks=100,main="",freq=FALSE,xlab="importances",ylab="density",ylim=c(0,1))
         curve(my.dsn(x,   xi=object$estimates[1], omega=object$estimates[2], lambda= object$estimates[3]),add=TRUE,col="red",lwd=2)
         abline(v=sn::qsn(0.95,  xi=object$estimates[1], omega=object$estimates[2], alpha= object$estimates[3]),col="red",lwd=2)
         abline(v=object$x[as.numeric(names(ww))],lwd=2,col="orange")
         abline(v= object$C,lwd=2,col="blue")
         abline(v= object$cc,lwd=2,col="purple")
         lines(object$x,object$fdr,lwd=3)
         abline(h=cutoff)
         legend("topright", c("fitted curve","95% quantile","cutoff","C","cc"),col=c("red","red", "orange","blue","purple"),lty=1,lwd=2)

     }
     if(debug.flag==2){
         cat(names(imp)[imp> object$x[as.numeric(names(ww))]],"\n")
     }
    
    a1<- match(names(imp)[imp> object$x[as.numeric(names(ww))]],names(imp))
    ppp<- 1-sn::psn(imp[a1], xi=object$estimates[1], omega=object$estimates[2], alpha= object$estimates[3])
    names(ppp) <-names(imp)[imp> object$x[as.numeric(names(ww))]]
     if (debug.flag==2){
         cat(length(ppp),"\n") 
     }

    ppp
}
