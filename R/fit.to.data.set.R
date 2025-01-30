#' fit.to.data.set
#'
#' This function fit a skew normal to a set of data
#' @param imp importances
#' @param df, contains x and y, midpoints and counts from a histogram of imp
#' @param debug.flag debug flag
#' @param plot.string, file name for a debugging plot
#' @param temp.dir, directory for debugging output
#' @param try.counter
#'        try.counter=1 my.dsn xi=  1  
#'        try.counter=2  xi=  mean(x)
#'        try.counter=3 start xi, omega, lambda from the parameters retuned by fitdistrplus::fitdist
#' @param return.all TRUE, return the full ouput of minpack.lm::nlsLM,
#'                   FALSE , return summary of parameters
#' @export
#' @importFrom graphics box legend lines hist
#' @importFrom  stats glm poisson
#' @importFrom grDevices dev.off png
#' @return If the skew-normal fitting routine is succesful, then the matrix of parmaters and standard errors is returned.
#'         -- othewise a "try-error" message is returned
#' @examples
#' \dontrun{
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
#' hist(temp, breaks = 200, freq = FALSE)
#' lines(df$x, df$y, type = "l", col = "green", lwd = 2, 
#'       xlim = c(0, max(df$x) + 0.5))
#' curve(sn::dsn(x, xi = fitted_parameters$Estimate[1], omega = fitted_parameters$Estimate[2], 
#'               alpha = fitted_parameters$Estimate[3]), add = TRUE, 
#'                 col = "purple", lwd = 3, xlim = c(0, 16))
#' curve(my.dsn(x, xi = fitted_parameters$Estimate[1], omega = fitted_parameters$Estimate[2],  
#'                 lambda = fitted_parameters$Estimate[3]), add = TRUE, 
#'                 col = "orange", lwd = 3)
#' }
#' 
#' \dontrun{
#' library(RFlocalfdr.data)
#' data(ch22)                                       
#' imp<-log(ch22$imp)                               
#' t2<-ch22$C                                       
#' temp<-imp[t2 > 30]   #                           
#' temp<-temp[temp != -Inf]                         
#' temp <- temp - min(temp) + .Machine$double.eps   
#' f_fit <- f.fit(temp)                             
#' y <- f_fit$zh$density                            
#' x <- f_fit$midpoints                             
#' df <- data.frame(x, y)                           
#' mm.df3 <- fit.to.data.set(df, temp)              
#' mm.df3
#' ##              Estimate Std..Error  t.value     Pr...t..
#' ## xi.xi        1.102303 0.03669284 30.04136 1.485263e-56
#' ## omega.omega  1.246756 0.04716184 26.43569 6.276349e-51
#' ## lambda.alpha 1.799169 0.17343872 10.37351 3.103195e-18
#' }

fit.to.data.set<-function (df, imp, debug.flag = 0, plot.string = "", temp.dir = NULL, 
    try.counter = 3, return.all = FALSE) 
{
    debug <- FALSE
    mm1.df <- NA
    class(mm1.df) <- "try-error"
    x <- df$x
    y <- df$y
    names(df) <- c("x", "y")
    if (inherits(mm1.df, "try-error") & try.counter == 1) {
        mm1.df <- try(minpack.lm::nlsLM(y ~ my.dsn(x, xi = xi, 
            omega = omega, lambda = lambda), start = list(xi = 1, 
            omega = 2, lambda = 1), data = df, control = minpack.lm::nls.lm.control(maxiter = 400)), 
            silent = TRUE)
        if (debug.flag > 0) {
            message(class(mm1.df), "class(mm1.df) -- try 1", "\n")
        }
    }
    if (inherits(mm1.df, "try-error") & (try.counter ==2)) {
        mm1.df <- try(minpack.lm::nlsLM(y ~ my.dsn(x, xi = xi, 
            omega = omega, lambda = lambda), start = list(xi = mean(x), 
            omega = 2, lambda = 1), data = df, control = minpack.lm::nls.lm.control(maxiter = 400)), 
            silent = TRUE)
        if (debug.flag > 0) {
            message(class(mm1.df), "class(mm1.df) -- try 2", "\n")
        }
    }
    if (inherits(mm1.df, "try-error") & try.counter == 3) {
        vip.sn.mle <- fitdistrplus::fitdist(imp, "sn", start = list(xi = mean(imp), 
                      omega = 1, alpha = 0), lower = c(-Inf, 0, -Inf), 
                      fix.arg = list(tau = 0), method = c("mle"))
        mm1.df <- try(minpack.lm::nlsLM(y ~ my.dsn(x, xi = xi, 
                   omega = omega, lambda = lambda), start = list(xi = vip.sn.mle$estimate[1], 
                   omega = vip.sn.mle$estimate[2], lambda = vip.sn.mle$estimate[3]), 
                   data = df, control = minpack.lm::nls.lm.control(maxiter = 400)), 
                   silent = TRUE)
        
        if (debug.flag > 0) {
            message(class(mm1.df), "class(mm1.df)-- try 3", "\n")
        }
    }

    if (debug.flag > 0 & !(inherits(mm1.df,"try-error"))) {
        png(paste(temp.dir, "/fit_to_data_set_", plot.string, 
            ".png", sep = ""))
        hist(imp, breaks = 200, freq = FALSE)
        lines(df$x, df$y, type = "l", col = "green", lwd = 2, 
            xlim = c(0, max(df$x) + 0.5))
        if (try.counter == 3) {
            curve(sn::dsn(x, xi = vip.sn.mle$estimate[1], omega = vip.sn.mle$estimate[2], 
                alpha = vip.sn.mle$estimate[3]), add = TRUE, 
                col = "purple", lwd = 3, xlim = c(0, 16))
            curve(my.dsn(x, xi = vip.sn.mle$estimate[1], omega = vip.sn.mle$estimate[2], 
                lambda = vip.sn.mle$estimate[3]), add = TRUE, 
                col = "purple", lwd = 3)
        }
        lines(x, predict(mm1.df), col = "red", lwd = 3)
        curve(my.dsn(x, xi = summary(mm1.df)$parameters[1, 1], 
            omega = summary(mm1.df)$parameters[2, 1], lambda = summary(mm1.df)$parameters[3, 
                1]), add = TRUE, col = "blue", lwd = 3)
        curve(sn::dsn(x, xi = summary(mm1.df)$parameters[1, 1], 
            omega = summary(mm1.df)$parameters[2, 1], alpha = summary(mm1.df)$parameters[3, 
                1]), add = TRUE, col = "blue", lwd = 3)
        abline(v = sn::qsn(0.5, xi = summary(mm1.df)$parameters[1, 
            1], omega = summary(mm1.df)$parameters[2, 1], alpha = summary(mm1.df)$parameters[3, 
            1]), col = "gray30")
        try(abline(v = sn::qsn(0.05, xi = summary(mm1.df)$parameters[1, 
            1], omega = summary(mm1.df)$parameters[2, 1], alpha = summary(mm1.df)$parameters[3, 
            1])), silent = TRUE)
        abline(v = sn::qsn(0.95, xi = summary(mm1.df)$parameters[1, 
            1], omega = summary(mm1.df)$parameters[2, 1], alpha = summary(mm1.df)$parameters[3, 
            1]), col = "gray30")
        message(sum(abs(df$y - predict(mm1.df))), " sum(abs(df$y-predict(mm1.df)))", 
            "\n")
        if (try.counter == 3) {
            legend("topright", c("spline fit", "fitted f0", "initial fitdist fit", 
                "Skew-normal at fitted values", "quantiles of skew normal"), 
                col = c("green", "red", "purple", "blue", "grey30"), 
                lty = 1, lwd = 4)
        }
        else {
            legend("topright", c("spline fit", "fitted f0", "Skew-normal at fitted values (my.dsn)",
                                 "Skew-normal at fitted values (dsn)",
                "quantiles of skew normal (dsn)"), col = c("green", 
                "red", "blue", "blue", "grey30"), lty = 1, lwd = c(4,4,4,4,1))
        }
        dev.off()
    }

    if (inherits(mm1.df,"try-error")){
        return.all <- TRUE
    }
    
    if (return.all == TRUE) {
        aa <- mm1.df
    }
    else {
         aa <- data.frame(summary(mm1.df)$parameters)
    }
    aa
}



