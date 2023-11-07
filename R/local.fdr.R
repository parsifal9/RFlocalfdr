#' local fdr 
#'
#' calculate the local 
#' @param f  object retuned by call to f.fit 
#' @param x  f_fit$midpoints    
#' @param FUN my.dsn
#' @param p0 estimated proportion of null importances
#' @param debug.flag  either 0 (no debugging information), 1 or 2
#' @param plot.string, file name for a debugging plot
#' @param temp.dir, directory for debugging output
#' @param ... arguments passed to FUN
#' @export
#' @return returns an estimate of the local false discovery rate.
#' @examples
#' data(imp20000)                                      
#' imp<-log(imp20000$importances)                               
#' t2<-imp20000$counts
#' temp<-imp[t2 > 1]   #see                          
#' temp<-temp[temp != -Inf]                         
#' temp <- temp - min(temp) + .Machine$double.eps   
#' f_fit <- f.fit(temp)                             
#' y <- f_fit$zh$density                            
#' x <- f_fit$midpoints                             
#' df <- data.frame(x, y)                           
#' fitted_parameters <- fit.to.data.set(df, temp, try.counter = 3)           
#' fitted_parameters
#' 
#' aa <- local.fdr(f_fit, df$x, FUN = my.dsn, xi = fitted_parameters$Estimate[1],
#'                 omega = fitted_parameters$Estimate[2], lambda = fitted_parameters$Estimate[3],
#'                 debug.flag = 0, plot.string = "initial")
#' 
#' plot(x,y,axes=FALSE,type="l",col="blue",main = "local fdr",
#'       xlab="importances",ylab="")                                                                                       
#' axis(2, pretty( c(0,max(y)+0.5*max(y)),10))                                                                            
#'                                                                                                                         
#' oldpar <- par(new = TRUE)
#' plot(x, aa, type="l",col="green",main = "",xlab="",ylab="",axes=FALSE)                                                 
#' abline(h = 0.2)                                                                                                        
#' axis(4, pretty( aa,10))                                                                                                
#'                                                                                                                         
#' axis(1,pretty(x,10))                                                                                                   
#' box() #- to make it look "as usual                                                                                     
#' legend("topright",c("density importances","local fdr"),col=c("blue","green"),lty=1)
#' par(oldpar)
#' 
#' \donttest{
#' library(RFlocalfdr.data)
#' data(ch22)                                                                                                             
#' imp<-log(ch22$imp)
#' t2<-ch22$C                                                                                                             
#' imp<-imp[t2 > 30]                                                                                                      
#' imp <- imp - min(imp) + .Machine$double.eps                                                                            
#' debug.flag <- 0                                                                                                        
#' f_fit <- f.fit(imp, debug.flag = debug.flag)                                                                           
#' y <- f_fit$zh$density                                                                                                  
#' x <- f_fit$midpoints                                                                                                   
#' df <- data.frame(x, y)                                                                                                 
#' initial.estimates <- fit.to.data.set.wrapper(df, imp, debug.flag = debug.flag,
#' return.all = FALSE)
#' 
#' aa <- local.fdr(f_fit, df$x, FUN = my.dsn, xi = initial.estimates$Estimate[1],
#'     omega = initial.estimates$Estimate[2], lambda = initial.estimates$Estimate[3],  debug.flag = 0,
#'                     plot.string = "initial")
#' 
#' plot(x,y,axes=FALSE,type="l",col="blue",main = "local fdr",                                                            
#'      xlab="importances",ylab="")                                                                                       
#' axis(2, pretty( c(0,max(y)+0.5*max(y)),10))                                                                            
#'                                                                                                                        
#' oldpar <- par(new = TRUE)
#' plot(x, aa, type="l",col="green",main = "",xlab="",ylab="",axes=FALSE)                                                 
#' abline(h = 0.2)                                                                                                        
#' axis(4, pretty( aa,10))                                                                                                
#'                                                                                                                        
#' axis(1,pretty(x,10))                                                                                                   
#' box() #- to make it look "as usual                                                                                     
#' legend("topright",c("density importances","local fdr"),col=c("blue","green"),lty=1)
#' par(oldpar)
#' }

local.fdr <-function(f,x,FUN=my.dsn, p0= 1, debug.flag=0,plot.string="",temp.dir=NULL, ...){
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


