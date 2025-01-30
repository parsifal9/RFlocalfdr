#' significant.genes
#'
#' This function sepects the significant "genes" and makes some plots
#' @param object object retruned by run.it.importance
#' @param imp importances
#' @param cutoff cutoff
#' @param do.plot do.plot either TRUE or FALSE (no plot)
#' @param debug.flag debug.flag  either 0 (no debugging information), 1 or 2 
#' @param use_95_q use the 0.95 q value
#' @keywords significant genes
#' @export
#' @return A list containg
#' - probabilities (from the fitted SN distribution) and names of the significant variables
#' - the estimated FDR 
#' @examples
#' \dontrun{
#' data(imp20000)
#' imp <- log(imp20000$importances)
#' t2  <- imp20000$counts
#' plot(density((imp)))
#' hist(imp,col=6,lwd=2,breaks=100,main="histogram of importances")
#' res.temp <- determine_cutoff(imp, t2 ,cutoff=c(0,1,2,3),plot=c(0,1,2,3),Q=0.75,try.counter=1)
#' plot(c(0,1,2,3),res.temp[,3])
#' imp<-imp[t2 > 1]
#' qq <- plotQ(imp,debug.flag = 0)                                                          
#' ppp<-run.it.importances(qq,imp,debug=0)                                                       
#' aa<-significant.genes(ppp,imp,cutoff=0.2,debug.flag=0,do.plot=2, use_95_q=TRUE)                           
#' length(aa$probabilities) #11#                                                          
#' names(aa$probabilities)
#' }
#'
#' \dontrun{
#' library(RFlocalfdr.data)
#' data(ch22)                                                                                 
#' ? ch22                                                                                     
#' plot(density(log(ch22$imp)))                                                               
#' t2 <-ch22$C                                                                                
#' imp<-log(ch22$imp)                                                                         
#' # Detemine a cutoff to get a unimodal density.                                              
#' # This may take several attempts. The default values of cutoff=c(0,1,4,10,15,20) will not find
#' # the minimum here.
#' #which occurs at 30
#' plot(c(25,30,35,40),res.temp[,3])                                                          
#' imp<-imp[t2 > 30]
#' qq <- plotQ(imp,debug.flag = 0)
#' ppp<-run.it.importances(qq,imp,debug=0)                                                       
#' aa<-significant.genes(ppp,imp,cutoff=0.2,debug.flag=0,do.plot=2)                           
#' length(aa$probabilities) # 6650                                                            
#' aa<-significant.genes(ppp,imp,cutoff=0.05,debug.flag=0,do.plot=2)                          
#' length(aa$probabilities) # 3653
#' }




significant.genes<-function (object, imp, cutoff = 0.2, use_95_q=TRUE, do.plot = TRUE, debug.flag = 0) 
{
    imp <- imp - min(imp) + .Machine$double.eps
    if (length(names(imp)) != length(imp)) {
        stop("imp should have names as this is what is returned\n")
    }
    start.x <- which.min(abs(object$x - mean(imp)))
    #plot(density(imp))
    #abline(v=mean(imp),col="red")
    #abline(v=object$x[start.x])
    
    #plot(object$x,object$fdr_0.95)
    #points(object$x[start.x:119],object$fdr_0.95[start.x:119],col="red" )
    #abline(h=cutoff)
    
    if (use_95_q){
        ww <- which.min(abs(object$fdr_0.95[start.x:119] - cutoff))
        if (do.plot == TRUE){do.plot=3}
    }
    else
    {
        ww <- which.min(abs(object$fdr_cc[start.x:119] - cutoff))
        if (do.plot == TRUE){do.plot=4}
    }
    
    #abline(v=object$x[as.numeric(names(ww))])
    
    num.sig.genes <- sum(imp > object$x[as.numeric(names(ww))])
    
    #do.plot == 1 not used
    #do.plot == 2 not used
    #do.plot == 3 fdr_0.95
    #do.plot == 2 fdr_cc
    oldpar <- par(no.readonly = TRUE )
    on.exit(par(oldpar))
    
    par(mar = c(4, 4, 4, 6))
    
    if (do.plot == 1) {
        aa <- hist(imp, col = 6, lwd = 2, breaks = 100, main = "", 
                   freq = FALSE, xlab = "importances", ylab = "density", 
                   axes = FALSE)
        axis(2, pretty(c(0, max(aa$density) + 0.5 * max(aa$density)), 
                       10))
        par(new = TRUE)
        plot(c(0, max(object$x)), c(0, 1), type = "n", axes = FALSE, 
             xlab = "", ylab = "")
        lines(object$x, object$fdr_0.95)
        abline(h = cutoff)
        abline(v = object$x[as.numeric(names(ww))])
        axis(4, pretty(c(0, 1), 10))
        mtext(side = 4, line = 3, "local fdr")
        axis(1, pretty(range(1:max(object$x)), 10))
        box()
    }
    if (debug.flag == 1) {
        message(sum(imp > object$x[as.numeric(names(ww))]), " sum(imp> x[as.numeric(names(ww))])","\n")
    }
    if (do.plot == 2) {
        aa <- hist(imp, col = 6, lwd = 2, breaks = 100, main = "", 
                   freq = FALSE, xlab = "importances", ylab = "density", 
                   axes = FALSE)
        x <- NULL; rm(x)
        curve(my.dsn(x, xi = object$estimates_C_0.95$Estimate[1], 
                     omega = object$estimates_C_0.95$Estimate[2], lambda = object$estimates_C_0.95$Estimate[3]), 
              add = TRUE, col = "red", lwd = 2)
        abline(v = sn::qsn(0.95, xi = object$estimates_C_0.95$Estimate[1], 
                           omega = object$estimates_C_0.95$Estimate[2], alpha = object$estimates_C_0.95$Estimate[3]), 
               col = "red", lwd = 2)
        abline(v = object$x[as.numeric(names(ww))], lwd = 2, 
               col = "orange")
        abline(v = object$C_0.95, lwd = 2, col = "blue")
        
        abline(v = object$cc, lwd = 2, col = "purple")
        axis(2, pretty(c(0, max(aa$density) + 0.5 * max(aa$density)), 
                       10))
        par(new = TRUE)
        plot(c(0, max(object$x)), c(0, 1), type = "n", axes = FALSE, 
             xlab = "", ylab = "")
        lines(object$x, object$fdr_0.95, lwd = 3)
        abline(h = cutoff)
        axis(4, pretty(c(0, 1), 10))
        mtext(side = 4, line = 3, "local fdr")
        axis(1, pretty(range(1:max(object$x)), 10))
        legend("topright", c("fitted curve", "95% quantile", 
                             "cutoff", "C", "cc", "fdr"), col = c("red", "red", 
                                                                  "orange", "blue", "purple", "black"), lty = 1, lwd = 2)
        box()
    }
    
    if (do.plot == 3) {
        aa <- hist(imp, col = 6, lwd = 2, breaks = 100, main = "", 
                   freq = FALSE, xlab = "importances", ylab = "density", 
                   axes = FALSE)
        x <- NULL; rm(x)
        curve(my.dsn(x, xi = object$estimates_C_0.95$Estimate[1], 
                     omega = object$estimates_C_0.95$Estimate[2], lambda = object$estimates_C_0.95$Estimate[3]), 
              add = TRUE, col = "red", lwd = 2)
        #    abline(v = sn::qsn(0.95, xi = object$estimates_C_0.95$Estimate[1], 
        #                       omega = object$estimates_C_0.95$Estimate[2], alpha = object$estimates_C_0.95$Estimate[3]), 
        #           col = "red", lwd = 2)
        abline(v = object$x[as.numeric(names(ww))], lwd = 2, 
               col = "orange")
        abline(v = object$C_0.95, lwd = 2, col = "blue")
        
        axis(2, pretty(c(0, max(aa$density) + 0.5 * max(aa$density)), 
                       10))
        par(new = TRUE)
        plot(c(0, max(object$x)), c(0, 1), type = "n", axes = FALSE, 
             xlab = "", ylab = "")
        lines(object$x, object$fdr_0.95, lwd = 3)
        abline(h = cutoff)
        axis(4, pretty(c(0, 1), 10))
        mtext(side = 4, line = 3, "local fdr")
        axis(1, pretty(range(1:max(object$x)), 10))
        legend("topright", c("fitted curve",  
                             "cutoff", "C",  "fdr"), col = c("red",  "orange", "blue", "black"), lty = 1, lwd = 2)
        box()
    }
    
    
    
    if (do.plot == 4) {
        aa <- hist(imp, col = 6, lwd = 2, breaks = 100, main = "", 
                   freq = FALSE, xlab = "importances", ylab = "density", 
                   axes = FALSE)
        x <- NULL; rm(x)
        curve(my.dsn(x, xi = object$estimates_cc$Estimate[1], 
                     omega = object$estimates_cc$Estimate[2], lambda = object$estimates_cc$Estimate[3]), 
              add = TRUE, col = "red", lwd = 2)
        abline(v = object$x[as.numeric(names(ww))], lwd = 2, 
               col = "orange")
        abline(v = object$cc, lwd = 2, col = "blue")
        
        axis(2, pretty(c(0, max(aa$density) + 0.5 * max(aa$density)), 
                       10))
        par(new = TRUE)
        plot(c(0, max(object$x)), c(0, 1), type = "n", axes = FALSE, 
             xlab = "", ylab = "")
        lines(object$x, object$fdr_cc, lwd = 3)
        abline(h = cutoff)
        axis(4, pretty(c(0, 1), 10))
        mtext(side = 4, line = 3, "local fdr")
        axis(1, pretty(range(1:max(object$x)), 10))
        legend("topright", c("fitted curve",  
                             "cutoff", "cc", "fdr"), col = c("red", "orange","blue", 
                                                             "black"), lty = 1, lwd = 2)
        box()
    }
    
    
    
    if (debug.flag == 2) {
        message(names(imp)[imp > object$x[as.numeric(names(ww))]], "\n")
    }
    a1 <- match(names(imp)[imp > object$x[as.numeric(names(ww))]],  names(imp))
    
    if (use_95_q){
        ppp <- 1 - sn::psn(imp[a1], xi = object$estimates_C_0.95$Estimat[1], 
                           omega = object$estimates_C_0.95$Estimat[2], alpha = object$estimates_C_0.95$Estimat[3])
        names(ppp) <- names(imp)[imp > object$x[as.numeric(names(ww))]]
        if (debug.flag == 2) {
            message(length(ppp), "\n")
        }
        cut <- 1 - sn::psn(object$x[as.numeric(names(ww))], xi = object$estimates_C_0.95$Estimat[1], 
                           omega = object$estimates_C_0.95$Estimat[2], alpha = object$estimates_C_0.95$Estimat[3])
        FDR <- cut * length(imp)/length(ppp)
        temp <- list(ppp, FDR)
        names(temp) <- c("probabilities", "FDR")
    }else {
        ppp <- 1 - sn::psn(imp[a1], xi = object$estimates_cc$Estimat[1], 
                           omega = object$estimates_cc$Estimat[2], alpha = object$estimates_cc$Estimat[3])
        names(ppp) <- names(imp)[imp > object$x[as.numeric(names(ww))]]
        if (debug.flag == 2) {
            message(length(ppp), "\n")
        }
        cut <- 1 - sn::psn(object$x[as.numeric(names(ww))], xi = object$estimates_cc$Estimat[1], 
                           omega = object$estimates_cc$Estimat[2], alpha = object$estimates_cc$Estimat[3])
        FDR <- cut * length(imp)/length(ppp)
        temp <- list(ppp, FDR)
        names(temp) <- c("probabilities", "FDR")
        
    }
    temp
}






