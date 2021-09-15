#' determine.C 
#'
#' This function allows you to express your love of cats.
#' @param f_fit object returned by f.fit
#' @param df data frame 
#' @param trace.plot Do you love cats? Defaults to TRUE.
#' @param starting_valu Do you love cats? Defaults to TRUE.
#' #' @keywords cats
#' @export
#' @examples
#' cat_function()
determine.C <-function (f_fit,df,trace.plot=FALSE,starting_value=1) 
{
    f <- f_fit$f.spline
    x<- f_fit$x
    y<-df$y
    qq <- rep(0, 119)
    if (starting_value==1){
        xi <- 1
        omega <- 2
        lambda <- 1
    }else{
        mm1.df = minpack.lm::nlsLM(y ~ my.dsn(x,xi=xi, omega=omega, lambda=lambda), 
                                   start = list(xi=1, omega=2, lambda= 1 ),
                                   data = df,control=minpack.lm::nls.lm.control(maxiter=400))
        
        xi <- summary(mm1.df)$parameters[1,1]
        omega <- summary(mm1.df)$parameters[2,1]
        lambda <-  summary(mm1.df)$parameters[3,1] 
    }
    
    for (ii in 15:119) {
        df2 <- df[1:ii, ]
        cat("dim(df2)", dim(df2), "\n")
        mm1.df2 = minpack.lm::nlsLM(y ~ my.dsn(x, xi = xi, omega = omega, 
            lambda = lambda), start = list(xi = xi, omega = omega, 
            lambda = lambda), data = df2, control = minpack.lm::nls.lm.control(maxiter = 400))
        xi <- summary(mm1.df2)$parameters[1, 1]
        omega <- summary(mm1.df2)$parameters[2, 1]
        lambda <- summary(mm1.df2)$parameters[3, 1]
        f0.1 <- predict(mm1.df2, newdata = df)
#        plot(x,f0.1)
#        f0.2 <- f0.1/(sum(diff(x)[1] * predict(mm1.df2, newdata = df))) #probably makes very li
#        lines(x,f0.2)
        ppp <- cumsum(f0.1) * diff(x)[1]

#        hist(ppp)
#        ppp<-sn::psn(x,  xi=summary(mm1.df2)$parameters[1,1], omega=summary(mm1.df2)$parameters[2,1], alpha= summary(mm1.df2)$parameters[3,1] )

        if (trace.plot==TRUE){
            plot(x,y,type="l",col="grey90",lwd=2,xlim=c(0,3))
            lines(df2$x,df2$y,col="green",lwd=2)
            lines(x,predict(mm1.df2,newdata=df),col="blue",lwd=3)  
            system("sleep 1")
        }

        p0 <- propTrueNullByLocalFDR(ppp)
        f0 <- (sum(f) * f0.1)/sum(f0.1)
        cat("p0 = ", p0, "\n")
        qq[ii] <- cumsum((-f_fit$counts * log(f0/(f))) - log(p0))[ii]
    }
    qq
}
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
