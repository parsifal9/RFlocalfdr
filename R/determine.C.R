#' determine.C 
#'
#' by assumption, there is a point $q$ such that to the left of $q$ $f_B \sim  f_0 (z)$. That is, there is a $q$
#' such that there are only null values to the left of $q$. We determine $q$ using a
#' change point method related to penalized model selection. See
#' Gauran, Iris Ivy M. and Park, Junyong and Lim, Johan and Park, DoHwan and Zylstra, John and Peterson,
#' Thomas and Kann, Maricel and Spouge, John L. "Empirical null estimation using zero-inflated discrete
#' mixture distributions and its application to protein domain data" Biometrics, 2018 74:2
#' @param f_fit object returned by f.fit
#' @param df data frame containg x and y
#' @param t1 initial estimates  xi.xi  omega.omega lambda. Probablgt returned by fit.to.data.set.wrapper
#' @param trace.plot -- produce a plot of each fit with a 1 second sleep. Can be watched as a movie.
#' @param starting_value -- needs discussion
#' @param start_at       -- needs discussion
#' @param debug             -- needs discussion
#' @keywords 
#' @export
#' @examples
#' data(ch22)                                                                                    
#' ? ch22                                                                                        
#' t2 <-ch22$C                                                                                   
#' imp<-log(ch22$imp)                                                                            
#' #Detemine a cutoff to get a unimodal density.                                                 
#' res.temp <- determine_cutoff(imp, t2 ,cutoff=c(25,30,35,40),plot=c(25,30,35,40),Q=0.75)       
#' plot(c(25,30,35,40),res.temp[,3])                                                             
#' imp<-imp[t2 > 30]
#'
#' f_fit<- f.fit(imp,debug.flag=debug.flag,temp.dir=temp.dir) #makes the plot histogram_of_variable_importances.png                              
#' y<-f_fit$zh$density                                                                                                                           
#' x<-f_fit$x                                                                                                                                    
#' plot(density(imp),main="histogram and fitted spline")                                                                                     
#' lines(x,y,col="red")                                                                                                                      
#' df<-data.frame(x,y)                                                                                                                           
#' initial.estimates <- fit.to.data.set.wrapper(df,imp,debug.flag=debug.flag,plot.string="initial",temp.dir=temp.dir,try.counter=try.counter)    
#' initial.estimates <- data.frame(summary(initial.estimates)$parameters)$Estimate                                                               
#'
#' qq<- determine.C(f_fit,df,initial.estimates,starting_value = 2,start_at=37,trace.plot = TRUE)    
#' cc<-x[which.min(qq)]                                                                             
#' plot(x,qq,main="determine cc")                                                                   
#' abline(v=cc)                                                                                     

determine.C<-function (f_fit, df, t1,trace.plot = FALSE, starting_value = 1,start_at=30,debug.flag=0) 
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
            cat("dim(df2)", dim(df2), "\n")
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
            plot(x, y, type = "l", col = "grey90", lwd = 2, xlim = c(0,  range(x)[[2]]+0.05))
            lines(df2$x, df2$y, col = "green", lwd = 2)
            lines(x, predict(mm1.df2, newdata = df), col = "blue",       lwd = 3)
            system("sleep 1")
        }
        p0 <- propTrueNullByLocalFDR(ppp)
        f0 <- (sum(f) * f0.1)/sum(f0.1)
        if(debug.flag > 0){
            cat("p0 = ", p0, "\n")
        }
        qq[ii] <- cumsum((-f_fit$counts * log(f0/(f))) - log(p0))[ii]
    }
    qq
}












## determine.C <-function (f_fit,df,trace.plot=FALSE,starting_value=1) 
## {
##     f <- f_fit$f.spline
##     x<- f_fit$x
##     y<-df$y
##     qq <- rep(0, 119)
##     if (starting_value==1){
##         xi <- 1
##         omega <- 2
##         lambda <- 1
##     }else{
##         mm1.df = minpack.lm::nlsLM(y ~ my.dsn(x,xi=xi, omega=omega, lambda=lambda), 
##                                    start = list(xi=1, omega=2, lambda= 1 ),
##                                    data = df,control=minpack.lm::nls.lm.control(maxiter=400))
        
##         xi <- summary(mm1.df)$parameters[1,1]
##         omega <- summary(mm1.df)$parameters[2,1]
##         lambda <-  summary(mm1.df)$parameters[3,1] 
##     }
    
##     for (ii in 15:119) {
##         df2 <- df[1:ii, ]
##         cat("dim(df2)", dim(df2), "\n")
##         mm1.df2 = minpack.lm::nlsLM(y ~ my.dsn(x, xi = xi, omega = omega, 
##             lambda = lambda), start = list(xi = xi, omega = omega, 
##             lambda = lambda), data = df2, control = minpack.lm::nls.lm.control(maxiter = 400))
##         xi <- summary(mm1.df2)$parameters[1, 1]
##         omega <- summary(mm1.df2)$parameters[2, 1]
##         lambda <- summary(mm1.df2)$parameters[3, 1]
##         f0.1 <- predict(mm1.df2, newdata = df)
## #        plot(x,f0.1)
## #        f0.2 <- f0.1/(sum(diff(x)[1] * predict(mm1.df2, newdata = df))) #probably makes very li
## #        lines(x,f0.2)
##         ppp <- cumsum(f0.1) * diff(x)[1]

## #        hist(ppp)
## #        ppp<-sn::psn(x,  xi=summary(mm1.df2)$parameters[1,1], omega=summary(mm1.df2)$parameters[2,1], alpha= summary(mm1.df2)$parameters[3,1] )

##         if (trace.plot==TRUE){
##             plot(x,y,type="l",col="grey90",lwd=2,xlim=c(0,3))
##             lines(df2$x,df2$y,col="green",lwd=2)
##             lines(x,predict(mm1.df2,newdata=df),col="blue",lwd=3)  
##             system("sleep 1")
##         }

##         p0 <- propTrueNullByLocalFDR(ppp)
##         f0 <- (sum(f) * f0.1)/sum(f0.1)
##         cat("p0 = ", p0, "\n")
##         qq[ii] <- cumsum((-f_fit$counts * log(f0/(f))) - log(p0))[ii]
##     }
##     qq
## }
## determine.C <-function(f_fit){
##     f<-f_fit$f.spline
##     qq<-rep(0,119)
##     xi<- 1
##     omega <- 2
##     lambda <- 1
##     for (ii in 15:119){
##         df2<-df[1:ii,]
##         cat("dim(df2)", dim(df2),"\n")
##         #    plot(df$x,df$y)
##         #    lines(df2$x,df2$y,col="red",lwd=2)
        
##         mm1.df2 = minpack.lm::nlsLM(y ~ my.dsn(x,xi=xi, omega=omega, lambda=lambda), 
##                                     start = list(xi=xi, omega=omega, lambda= lambda ),
##                                     data = df2,control=minpack.lm::nls.lm.control(maxiter=400))
##         xi <- summary(mm1.df2)$parameters[1,1]     
##         omega <- summary(mm1.df2)$parameters[2,1]
##         alpha <-  summary(mm1.df2)$parameters[3,1]
        
##         #   hist(imp1,breaks=100,freq=FALSE)
##         f0<- predict(mm1.df2,newdata=df)
##         f0 <-f0/(sum(diff(x)[1]*predict(mm1.df2,newdata=df)) )
##         #   lines(x,f0)    
##         ppp<- cumsum(f0)*diff(x)[1]
##         p0<- propTrueNullByLocalFDR(ppp)    
##         f0 <- (sum(f) * f0)/sum(f0)
##         cat("p0 = ",p0,"\n")
##         qq[ii]<- cumsum((-f_fit$counts*log(f0/(f))) - log(p0))[ii]
##     }
##     qq
## }
