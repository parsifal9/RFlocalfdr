#' determine.C 
#'
#' by assumption, there is a point q such that to the left of q, f_B sim  f_0 (z). That is, there is a q
#' such that there are only null values to the left of q. We determine q using a
#' change point method related to penalized model selection. See
#' Gauran, Iris Ivy M. and Park, Junyong and Lim, Johan and Park, DoHwan and Zylstra, John and Peterson,
#' Thomas and Kann, Maricel and Spouge, John L. "Empirical null estimation using zero-inflated discrete
#' mixture distributions and its application to protein domain data" Biometrics, 2018 74:2
#' @param f_fit object returned by f.fit
#' @param df data frame containing x and y
#' @param t1 initial estimates of xi, omega, and  lambda. Generally returned by fit.to.data.set.wrapper
#' @param trace.plot -- produce a plot of each fit with a 1 second sleep. Can be watched as a movie.
#' @param start_at       --  x <- f_fit$midpoints  is of length 119 (quite arbitrary). We use the first start_at  
#'                          values of x to fit the skew-normal distribution. 
#' @param debug.flag     -- debugging level. If debug.flag >0 then some output is printed to the screen. 
#' @importFrom graphics abline axis box curve legend lines mtext par
#' @importFrom  stats density predict quantile
#' @export
#' @return -- a vector of numbers of length equal to the rows in df (119 in this case). Say that this is qq.
#'           We determine the minimum value of qq. This is the value "C" such that
#'           -- to the right of C, our data is generated from the NULL distribution
#'           -- to the left of C, we have a mixture of the NULL and non-NULL distribution
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
#' initial.estimates <- fit.to.data.set.wrapper(df, temp, try.counter = 3,return.all=FALSE)           
#' initial.estimates<-  initial.estimates$Estimate
#' 
#' qq<- determine.C(f_fit,df,initial.estimates,start_at=37,trace.plot = FALSE)    
#' cc<-x[which.min(qq)]                                                                             
#' plot(x,qq,main="determine cc")                                                                   
#' abline(v=cc)
#' # unfortunately the minima does not appear reasonable. In this case it is advisable to use the
#' # 95th quantile
#' 
#' \donttest{
#' #needs the  chromosome 22 data in  RFlocalfdr.data. Also has a long runtime.
#' library(RFlocalfdr.data)
#' data(ch22)                                                                                    
#' ?ch22                                                                                        
#' t2 <-ch22$C                                                                                   
#' imp<-log(ch22$imp)                                                                            
#' #Detemine a cutoff to get a unimodal density.                                                 
#' res.temp <- determine_cutoff(imp, t2 ,cutoff=c(25,30,35,40),plot=c(25,30,35,40),Q=0.75)       
#' plot(c(25,30,35,40),res.temp[,3])                                                             
#' imp<-imp[t2 > 30]
#' debug.flag <- 0
#' f_fit<- f.fit(imp,debug.flag=debug.flag,temp.dir=temp.dir)
#' #makes the plot histogram_of_variable_importances.png                              
#' y<-f_fit$zh$density                                                                                                                           
#' x<-f_fit$midpoints                                                                                                                                    
#' plot(density(imp),main="histogram and fitted spline")                                                                                     
#' lines(x,y,col="red")                                                                                                                      
#' df<-data.frame(x,y)                                                                                                                           
#' initial.estimates <- fit.to.data.set.wrapper(df,imp,debug.flag=debug.flag,plot.string="initial",
#'                                               temp.dir=temp.dir,try.counter=3)    
#' initial.estimates <- data.frame(summary(initial.estimates)$parameters)$Estimate                                                               
#' # 1.102303 1.246756 1.799169
#' qq<- determine.C(f_fit,df,initial.estimates,start_at=37,trace.plot = TRUE)    
#' cc<-x[which.min(qq)]                                                                             
#' plot(x,qq,main="determine cc")                                                                   
#' abline(v=cc)
#' }

determine.C<-function (f_fit, df, t1, trace.plot = FALSE ,start_at=30, debug.flag=0) 
{
    f <- f_fit$f.spline
    x <- df$x
    y <- df$y
    qq <- rep(NA, 119)
    ## if (starting_value == 1) {
    ##     xi <- 1
    ##     omega <- 2
    ##     lambda <- 1
    ## }
    ## else {
    ##     mm1.df = minpack.lm::nlsLM(y ~ my.dsn(x, xi = xi, omega = omega, 
    ##         lambda = lambda), start = list(xi = mean(x), omega = 2, 
    ##         lambda = 1), data = df, control = minpack.lm::nls.lm.control(maxiter = 400))
    ##     xi <- summary(mm1.df)$parameters[1, 1]
    ##     omega <- summary(mm1.df)$parameters[2, 1]
    ##     lambda <- summary(mm1.df)$parameters[3, 1]
    ## }

    for (ii in start_at:119) {
        df2 <- df[1:ii, ]
        if(debug.flag > 0){
            message("dim(df2)", dim(df2), "\n")
            }
        mm1.df2 = minpack.lm::nlsLM(y ~ my.dsn(x, xi = xi, omega = omega, 
            lambda = lambda), start = list(xi = t1[1], omega = t1[2], 
            lambda = t1[3]), data = df2, control = minpack.lm::nls.lm.control(maxiter = 400))
        xi <- summary(mm1.df2)$parameters[1, 1]
        omega <- summary(mm1.df2)$parameters[2, 1]
        lambda <- summary(mm1.df2)$parameters[3, 1]
        f0.1 <- predict(mm1.df2, newdata = df)
        f0.1 <- f0.1 + .Machine$double.eps
        ppp <- cumsum(f0.1) * diff(x)[1]
        if (trace.plot == TRUE) {
             plot(x, y, type = "l", col = "grey90", lwd = 2 , xlim = c(range(x)[[1]]+0.05,  range(x)[[2]]+0.05))
            lines(df2$x, df2$y, col = "green", lwd = 2)
            lines(x, predict(mm1.df2, newdata = df), col = "blue",       lwd = 3)
            system("sleep 1")
        }
        p0 <- propTrueNullByLocalFDR(ppp)
        f0 <- (sum(f) * f0.1)/sum(f0.1)
        if(debug.flag > 0){
            message("p0 = ", p0, "\n")
        }
        qq[ii] <- cumsum((-f_fit$counts * log(f0/(f))) - log(p0))[ii]
    }
    qq
}










