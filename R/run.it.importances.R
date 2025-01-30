#' run.it.importances
#'
#' @param qq object retuned by plotQ
#' @param imp "reduction in impurity" importances from a random forest model
#' @param debug.flag  either 0 (no debugging information), 1 or 2
#' @param temp.dir if debug flag is >0 then information is written to temp.dir
#' @keywords variable importance 
#' @export
#' @return return a list contining
#' - "C_0.95" estimate of the cutoff "C" such that there are only null values to the left of C.
#'            Based on the 95th quantile of the density
#' - "cc"    estimate of the cutoff "C" based on the procedure of Gauran, et al., 2918, Biometrics,
#' - "estimates_C_0.95" estimate of the parameters of the SN using the data up to the C estimate
#' - "estimates_cc"     estimate of the parameters of the SN using the data up to the C estimate
#' - "fdr_0.95"       estimate of the fdr curve using the SN from "estimates_C_0.95"
#' - "fdr_cc"         estimate of the fdr curve using the SN from "estimates_cc"
#' - "x"        the x values from plotQ
#' - "temp.dir" the temp directory for dubuggin
#' - "p0"     the estmate of the proportion of null values (can be 1)
#' @examples
#' \dontrun{
#' data(imp20000)
#' imp <- log(imp20000$importances)
#' t2 <- imp20000$counts
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
#' #document how the data set is created                                                      
#' plot(density(log(ch22$imp)))                                                               
#' t2 <-ch22$C                                                                                
#' imp<-log(ch22$imp)                                                                         
#' #Detemine a cutoff to get a unimodal density.                                              
#' res.temp <- determine_cutoff(imp, t2 ,cutoff=c(1,10,20,30),plot=c(1,10,20,30),Q=0.75)      
#' plot(c(1,2,3,4),res.temp[,3])                                                              
#'                                                                                            
#' res.temp <- determine_cutoff(imp, t2 ,cutoff=c(25,30,35,40),plot=c(25,30,35,40),Q=0.75)    
#' plot(c(25,30,35,40),res.temp[,3])                                                          
#' imp<-imp[t2 > 30]
#' qq <- plotQ(imp,debug.flag = 0)                                                          
#' ppp<-run.it.importances(qq,imp,debug=0)                                                       
#' aa<-significant.genes(ppp,imp,cutoff=0.2,debug.flag=0,do.plot=2)                           
#' length(aa$probabilities) # 6650                                                            
#' aa<-significant.genes(ppp,imp,cutoff=0.05,debug.flag=0,do.plot=2)                          
#' length(aa$probabilities) # 3653
#' }


run.it.importances <- function(qq, imp, debug.flag = 0, temp.dir = NULL) {
    x <- qq$df$x
    y <- qq$df$y
    df <- qq$df
    final.estimates_C_0.95 <- qq$final.estimates_C_0.95
    final.estimates_cc <- qq$final.estimates_cc
    temp.dir <- temp.dir
    C_0.95 <- qq$C_0.95
    cc <- qq$cc
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    
    imp <- imp - min(imp) + .Machine$double.eps
    if (!is.na(cc)) {
        df3 <- data.frame(x[x < cc], y[x < cc])
        names(df3) <- c("x", "y")
    }
    if (debug.flag > 0) {
        if (length(temp.dir) == 0) {
            temp.dir <- tempdir()
        }
        fileConn <- file(paste(temp.dir, "/output_from_run_it_importances.txt", 
            sep = ""), open = "wt")
        writeLines(c("Hello", "World"), fileConn)
    }
    df2 <- data.frame(x[x < C_0.95], y[x < C_0.95])
    f_fit <- qq$f_fit
    ww <- qq$ww
    if (debug.flag > 0) {
        png(paste(temp.dir, "/compare_C_and_cc_and_the_resulting_fits.png", 
            sep = ""))
        hist(imp, col = 6, lwd = 2, breaks = 100, main = "compare C and cc and the resulting fits")
        abline(v = C_0.95, col = "red", lwd = 2)
        if (!is.na(cc)) {
            abline(v = cc, col = "purple", lwd = 2)
        }
        legend("topright", c("C", "cc"), col = c("red", "purple"), 
            lty = 1)
        dev.off()
        png(paste(temp.dir, "/compare_C_and_cc_and_the_resulting_fits_2.png", 
            sep = ""))
        hist(imp, breaks = 200, freq = FALSE)
        lines(x, y, type = "l", col = "grey90", lwd = 2, xlim = c(0, 
            12))
        lines(df2$x, df2$y, col = "green", lwd = 2)
        abline(v = C_0.95, col = "green")
        curve(sn::dsn(x, xi = final.estimates_C_0.95$Estimate[1], 
            omega = final.estimates_C_0.95$Estimate[2], alpha = final.estimates_C_0.95$Estimate[3]), 
            add = TRUE, col = "green", lwd = 3)
        lines(x, my.dsn(x, xi = final.estimates_C_0.95$Estimate[1], 
            omega = final.estimates_C_0.95$Estimate[2], lambda = final.estimates_C_0.95$Estimate[3]), 
            col = "green", lwd = 3)
        if (!is.na(cc)) {
            lines(df3$x, df3$y, col = "blue", lwd = 2)
            abline(v = cc, col = "blue")
            if (!inherits(final.estimates_cc,"character")) {
                curve(sn::dsn(x, xi = final.estimates_cc$Estimate[1], 
                  omega = final.estimates_cc$Estimate[2], alpha = final.estimates_cc$Estimate[3]), 
                  add = TRUE, col = "blue", lwd = 3)
            }
        }
        legend("topright", c("C", "cc"), col = c("green", "blue"), 
            lty = 1)
        dev.off()
        png(paste(temp.dir, "/fit.to.data.set_df2.png", sep = ""))
        plot(x, y, type = "l", col = "grey90", lwd = 2)
        lines(df2$x, df2$y, col = "green", lwd = 2)
        lines(x, my.dsn(x, xi = final.estimates_C_0.95$Estimate[1], 
            omega = final.estimates_C_0.95$Estimate[2], lambda = final.estimates_C_0.95$Estimate[3]))
        dev.off()
    }
    ppp <- sn::psn(imp, xi = final.estimates_C_0.95$Estimate[1], 
        omega = final.estimates_C_0.95$Estimate[2], alpha = final.estimates_C_0.95$Estimate[3])
    p0 <- propTrueNullByLocalFDR(ppp)
    if (debug.flag == 1) {
        writeLines(paste(p0, "p0"), fileConn)
    }
    aa_C_0.95 <- local.fdr(f_fit, df$x, FUN = my.dsn, xi = final.estimates_C_0.95$Estimate[1], 
        omega = final.estimates_C_0.95$Estimate[2], lambda = final.estimates_C_0.95$Estimate[3], 
        p0 = p0, debug.flag = debug.flag, plot.string = "final", 
        temp.dir = temp.dir)
    aa_cc <- try(local.fdr(f_fit, df$x, FUN = my.dsn, xi = final.estimates_cc$Estimate[1], 
        omega = final.estimates_cc$Estimate[2], lambda = final.estimates_cc$Estimate[3], 
        p0 = p0, debug.flag = debug.flag , plot.string = "final", 
        temp.dir = temp.dir), silent = TRUE)
    if (debug.flag > 0) {
        png(paste(temp.dir, "/run.it.importances_final_plot.png", 
            sep = ""))
        par(mar = c(4, 4, 4, 6))
        aa <- hist(imp, col = 6, lwd = 2, breaks = 100, main = "", 
            freq = FALSE, xlab = "importances", ylab = "density", 
            axes = FALSE)
        abline(v = c(1:100))
        curve(my.dsn(x, xi = final.estimates_C_0.95$Estimate[1], 
            omega = final.estimates_C_0.95$Estimate[2], lambda = final.estimates_C_0.95$Estimate[3]), 
            add = TRUE, col = "red", lwd = 2)
        abline(v = sn::qsn(0.95, xi = final.estimates_C_0.95$Estimate[1], 
            omega = final.estimates_C_0.95$Estimate[2], alpha = final.estimates_C_0.95$Estimate[3]), 
            col = "red", lwd = 2)
        if (!inherits(final.estimates_cc,"character")) {
            curve(my.dsn(x, xi = final.estimates_cc$Estimate[1], 
                omega = final.estimates_cc$Estimate[2], lambda = final.estimates_cc$Estimate[3]), 
                add = TRUE, col = "lawngreen", lwd = 2)
            abline(v = sn::qsn(0.95, xi = final.estimates_cc$Estimate[1], 
                omega = final.estimates_cc$Estimate[2], alpha = final.estimates_cc$Estimate[3]), 
                col = "lawngreen", lwd = 2)
        }
        abline(v = x[as.numeric(names(ww))], lwd = 2, col = "orange")
        abline(v = C_0.95, lwd = 2, col = "blue")
        abline(v = cc, lwd = 2, col = "purple")
        axis(2, pretty(c(0, max(aa$density) + 0.5 * max(aa$density)), 
            10))
        par(new = TRUE)
        plot(c(0, max(x)), c(0, 1), type = "n", axes = FALSE, 
            xlab = "", ylab = "")
        lines(x, aa_C_0.95, lwd = 3)
        abline(h = 0.2, col = "orange")
        axis(4, pretty(c(0, 1), 10))
        mtext(side = 4, line = 3, "local fdr")
        axis(1, pretty(range(1:max(x)), 10))
        legend("topright", c("fitted curve", "95% quantile", 
            "cutoff", "C", "cc", "fdr"), col = c("red", "red", 
            "orange", "blue", "purple", "black"), lty = 1, lwd = 2)
        box()
        dev.off()
    }
    if (debug.flag > 1) {
        message(names(imp)[imp > x[as.numeric(names(ww))]], "\n\n\n\n")
    }
    if (debug.flag > 0) {
        writeLines(c("GoodBye World"), fileConn)
        close(fileConn)
    }
    temp <- list(aa_C_0.95, aa_cc, df$x, final.estimates_C_0.95, 
        final.estimates_cc, temp.dir, C_0.95, cc, p0)
    names(temp) <- c("fdr_0.95", "fdr_cc", "x", "estimates_C_0.95,", 
        "estimates_cc", "temp.dir", "C_0.95", "cc", "p0")
    temp
}



