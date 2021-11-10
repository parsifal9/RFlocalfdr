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
        old.par <- par(mar = c(2,2,2,2))
         ## do plotting stuff with new settings
         par(mar = c(3, 3, 3, 6)) # Set the margin on all sides to 2
         aa<- hist(imp, col = 6, lwd = 2, breaks = 100, main = "", 
                   freq = FALSE, xlab = "importances", ylab = "density", 
                   axes=FALSE)  #ylim = c(0, max(aa$density)+0.5*max(aa$density))
         axis(2, pretty( c(0,max(aa$density)+0.5*max(aa$density)),10))
         par(new=TRUE)
         plot(c(0,max(object$x)),c(0,1),type="n",axes=FALSE, xlab = "", ylab = "")
         lines(object$x, object$fdr)
         abline(h = cutoff)
         abline(v = object$x[as.numeric(names(ww))])
         axis(4, pretty(c(0,1),10))
         mtext(side=4,line=3,  "local fdr")
         axis(1,pretty(range(1:max(object$x)),10))
         box() #- to make it look "as usual
         par(old.par)
     }
    
    if(debug.flag==1){
        cat(sum(imp> object$x[as.numeric(names(ww))]),"sum(imp> x[as.numeric(names(ww))])","\n")
    }
     if (do.plot==2){
         old.par <- par(mar = c(2,2,2,2))
         ## do plotting stuff with new settings
         par(mar = c(3, 3, 3, 6)) # Set the margin on all sides to 2
         aa<- hist(imp, col = 6, lwd = 2, breaks = 100, main = "", 
                   freq = FALSE, xlab = "importances", ylab = "density",axes=FALSE)
         curve(my.dsn(x, xi = object$estimates[1], omega = object$estimates[2], 
                      lambda = object$estimates[3]), add = TRUE, col = "red", 
               lwd = 2)
         abline(v = sn::qsn(0.95, xi = object$estimates[1], omega = object$estimates[2], 
                            alpha = object$estimates[3]), col = "red", lwd = 2)
         abline(v = object$x[as.numeric(names(ww))], lwd = 2, 
                col = "orange")
         abline(v = object$C, lwd = 2, col = "blue")
         abline(v = object$cc, lwd = 2, col = "purple")
         
         axis(2, pretty( c(0,max(aa$density)+0.5*max(aa$density)),10))
         par(new=TRUE)
         plot(c(0,max(object$x)),c(0,1),type="n",axes=FALSE, xlab = "", ylab = "")
         lines(object$x, object$fdr, lwd = 3)
         abline(h = cutoff)
         abline(v = object$x[as.numeric(names(ww))])
         axis(4, pretty(c(0,1),10))
         mtext(side=4,line=3,  "local fdr")
         axis(1,pretty(range(1:max(object$x)),10))
         legend("topright", c("fitted curve", "95% quantile", 
                              "cutoff", "C", "cc"), col = c("red", "red", "orange", 
                                                            "blue", "purple"), lty = 1, lwd = 2)
         
         box() #- to make it look "as usual
         par(old.par)
         
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
