#' plotQ
#'
#' produces a plot showing the q values 
#' * q_95, the 95th quantile of the data
#' * q using the penalized selection method of Gauran et.al 2018
#' 
#' We estiamte a value "q" such that:
#' to the left of "q", the density is composed solely of NULL importance values
#' to the right of "q" we have  a density that is a mixture of null and non-null  importance values.
#' The method of Gauran et.al 2018 may not work in cases where the data distribution is not well modelled by a skew-normal.
#' The q_95 value can be uses as a workaround in these case.
#' In many cases they will be very similar
#' 
#' @param imp "reduction in impurity" importances from a random forest model
#' @param debug.flag  either 0 (no debugging information), 1 or 2
#' @param temp.dir if debug flag is >0 then information is written to temp.dir
#' @param try.counter where to explain this?
#' @keywords variable importance
#' @return
#' - df, contains x and y, midpoints and counts from a histogram of imp
#' - final.estimates_C_0.95, the output from the fitting routine nlsLM in minpack.lm, where df has been truncated
#'   at the value C_0.95 (the 0.95 quantile of the skew-Normal distribution fitted to the imp histogram
#' - final.estimates_cc, as for final.estimates_C_0.95 but with cc determined by the procedure of Gauran et.al 2018
#' - temp.dir,  the directory where debugging information may be written
#' - C_0.95, the 0.95 quantile of the skew-Normal distribution fitted to the imp histogram
#' - cc,  determined by the procedure of Gauran et.al 2018
#' - fileConn, a file connectin for writing debugging information
#' - f_fit, a spline fit to the histogram
#' - ww the minimum value of the local fdr
#' @export
#' @examples
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
#'
#' \donttest{
#' library(RFlocalfdr.data)
#' data(ch22)                                                                                 
#' ?ch22                                                                                     
#' #document how the data set is created                                                      
#' plot(density(log(ch22$imp)))                                                               
#' t2 <-ch22$C                                                                                
#' imp<-log(ch22$imp)                                                                         
#' #Detemine a cutoff to get a unimodal density.
#' # This was calculated previously. See determine_cutoff
#' imp<-imp[t2 > 30]
#' qq <- plotQ(imp,debug.flag = 0)
#' 
#' data(smoking)
#' ?smoking 
#' y<-smoking$y
#' smoking_data<-smoking$rma
#' y.numeric <-ifelse((y=="never-smoked"),0,1)
#' 
#' library(ranger)
#' rf1 <-ranger::ranger(y=y.numeric ,x=smoking_data,importance="impurity",seed=123, num.trees = 10000,
#'              classification=TRUE)
#' t2 <-count_variables(rf1)
#' imp<-log(rf1$variable.importance)
#' plot(density(imp),xlab="log importances",main="")
#' cutoffs <- c(2,3,4,5)
#' res.con<- determine_cutoff(imp,t2,cutoff=cutoffs,plot=c(2,3,4,5))
#' 
#' plot(cutoffs,res.con[,3],pch=15,col="red",cex=1.5,ylab="max(abs(y - t1))")
#' cutoffs[which.min(res.con[,3])]
#' 
#' temp<-imp[t2 > 3]
#' temp <- temp - min(temp) + .Machine$double.eps
#' qq <- plotQ(temp)
#' ppp<-run.it.importances(qq,temp,debug.flag = 0)
#' aa<-significant.genes(ppp,temp,cutoff=0.05,debug.flag=0,do.plot=TRUE,use_95_q=TRUE)
#' length(aa$probabilities) # 17
#' }

plotQ <- function (imp, debug.flag = 0, temp.dir = NULL, try.counter = 3){
    fileConn <- NULL
    ww <- NULL
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    if (debug.flag > 0) {
        if (length(temp.dir) == 0) {
            temp.dir <- tempdir()
        }
        fileConn <- file(paste(temp.dir, "/output_from_PlotQ.txt", sep = ""), open = "wt")
        writeLines(c("Hello", "World"), fileConn)
    }
    imp <- imp - min(imp) + .Machine$double.eps
    if (debug.flag > 0) {
        writeLines(paste(length(imp), "length(imp)"), fileConn)
        writeLines(paste(unlist(range(imp)), "range(imp)"), fileConn)
    }
    f_fit <- f.fit(imp, debug.flag = debug.flag, temp.dir = temp.dir)
    y <- f_fit$zh$density
    x <- f_fit$midpoints
    if (debug.flag > 0) {
        png(paste(temp.dir, "/density_importances.png", sep = ""))
        plot(density(imp), main = "histogram and fitted spline")
        lines(x, y, col = "red")
        dev.off()
    }
    df <- data.frame(x, y)
    initial.estimates <- fit.to.data.set.wrapper(df, imp, debug.flag = debug.flag,
        plot.string = "initial", temp.dir = temp.dir, try.counter = try.counter)
    if (inherits(initial.estimates, "try-error")) {
        stop("fit.to.data.set.wrapper has returned a try-error at the initial estimate")
    } else{
        initial.estimates <- data.frame(summary(initial.estimates)$parameters)$Estimate
    }
    if (debug.flag > 0) {
        writeLines(paste("initial estimates", initial.estimates[1],
          initial.estimates[2], initial.estimates[3]), fileConn)
        writeLines("we calcualte the fdr using the initial estimates",
          fileConn)
        aa <- local.fdr(f_fit, df$x, FUN = my.dsn, xi = initial.estimates[1],
                        omega = initial.estimates[2], lambda = initial.estimates[3],
            debug.flag = 0, plot.string = "initial", temp.dir = temp.dir)
        png(paste(temp.dir, "/fdr_using_initial_estiamtes.png",
           sep = ""))
        plot(x, aa, main = "fdr using initial estimates")
        abline(h = 0.2)
        ww <- which.min(abs(aa[50:119] - 0.2))
        sum(imp > x[as.numeric(names(ww))])
        tt <- sum(imp > x[as.numeric(names(ww))])
        abline(v = x[as.numeric(names(ww))])
        dev.off()
        writeLines(paste("sum(imp> x[as.numeric(names(ww))])",
           tt), fileConn)
    }

    C_0.95 <- sn::qsn(0.95, xi = initial.estimates[1], omega = initial.estimates[2],
       alpha = initial.estimates[3])
    if (debug.flag == 1) {
        writeLines(paste("calculating C_0.95", C_0.95), fileConn)
    }
    if (debug.flag > 0) {
        message("calculating cc", "\n")
    }
    qq <- try(determine.C(f_fit, df, initial.estimates, start_at = 37), silent = TRUE)
    if (debug.flag > 0) {
        writeLines(paste(class(qq), "class(determine.C)"), fileConn)
    }
    cc <- final.estimates_cc <- aa_cc <- NA
    if (inherits(qq, "numeric")) {
        cc <- x[which.min(qq)]
        if (debug.flag > 0) {
            writeLines(paste("cc= ", cc), fileConn)
            png(paste(temp.dir, "/determine_cc.png", sep = ""))
            plot(x, qq, main = "determine cc")
            abline(v = cc)
            dev.off()
        }
    }
    df2 <- data.frame(x[x < C_0.95], y[x < C_0.95])
    names(df2) <- c("x", "y")
    final.estimates_C_0.95 <- fit.to.data.set.wrapper(df2, imp,
           debug.flag = debug.flag, plot.string = "final", temp.dir = temp.dir, try.counter = try.counter)
    if (inherits(final.estimates_C_0.95, "try-error")) {
        if (debug.flag > 0){
            print("fit.to.data.set.wrapper has returned a try-error at the final 0.95 estimates")
        }
    }
 
    final.estimates_C_0.95 <- data.frame(summary(final.estimates_C_0.95)$parameters)
    if (!is.na(cc)) {
        df3 <- data.frame(x[x < cc], y[x < cc])
        names(df3) <- c("x", "y")
        final.estimates_cc <- fit.to.data.set.wrapper(df3, imp, 
            debug.flag = debug.flag, plot.string = "cc", temp.dir = temp.dir)
        if  (inherits(final.estimates_cc, "try-error")) {
             if (debug.flag > 0){
                 print("fit.to.data.set.wrapper has returned a try-error at the final cc estimates")
             }
        } else {
            final.estimates_cc <- data.frame(summary(final.estimates_cc)$parameters)
        }
    }
    if (1) {
        aa <- hist(imp, col = "grey", lwd = 2, breaks = 100, 
            main = "", freq = FALSE, xlab = "importances", ylab = "density", 
            axes = FALSE)
        abline(v = C_0.95, col = "red", lwd = 2)
        if (!is.na(cc)) {
            abline(v = cc, col = "purple", lwd = 2)
        }
        legend("topright", c("q_95", "q"), col = c("red", "purple"), 
            lty = 1)
        axis(2, pretty(c(0, max(aa$density) + 0.5 * max(aa$density)), 
            10))
        lines(x, y, type = "l", col = "grey90", lwd = 2, xlim = c(0, 
            12))
        lines(df2$x, df2$y, col = "red", lwd = 2)
        if (!is.na(cc)) {
            lines(df3$x, df3$y, col = "purple", lwd = 2)
        }
        if (inherits(qq, "numeric")) {
            par(new = TRUE)
            plot(x, qq, type = "l", lwd = 2, axes = FALSE, xlab = "", 
                ylab = "", col = "purple")
            axis(4, pretty(c(min(qq, na.rm = TRUE), max(qq, na.rm = TRUE) + 
                0.5 * max(qq, na.rm = TRUE)), 10))
        }
        box()
        temp <- list(df, final.estimates_C_0.95, final.estimates_cc, 
            temp.dir, C_0.95, cc, f_fit, ww)
        names(temp) <- c("df", "final.estimates_C_0.95,", "final.estimates_cc", 
            "temp.dir", "C_0.95", "cc", "f_fit", "ww")
    }
    if (debug.flag > 0) {
        writeLines(c("GoodBye World"), fileConn)
        close(fileConn)
    }
    temp
}

