#' determine.C 
#'
#' This function allows you to express your love of cats.
#' @param f_fit Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
determine.C <-function(f_fit){
    f<-f_fit$f.spline
    qq<-rep(0,119)
    xi<- 1
    omega <- 2
    lambda <- 1
    for (ii in 15:119){
        df2<-df[1:ii,]
        cat("dim(df2)", dim(df2),"\n")
        #    plot(df$x,df$y)
        #    lines(df2$x,df2$y,col="red",lwd=2)
        
        mm1.df2 = minpack.lm::nlsLM(y ~ my.dsn(x,xi=xi, omega=omega, lambda=lambda), 
                                    start = list(xi=xi, omega=omega, lambda= lambda ),
                                    data = df2,control=minpack.lm::nls.lm.control(maxiter=400))
        xi <- summary(mm1.df2)$parameters[1,1]     
        omega <- summary(mm1.df2)$parameters[2,1]
        alpha <-  summary(mm1.df2)$parameters[3,1]
        
        #   hist(imp1,breaks=100,freq=FALSE)
        f0<- predict(mm1.df2,newdata=df)
        f0 <-f0/(sum(diff(x)[1]*predict(mm1.df2,newdata=df)) )
        #   lines(x,f0)    
        ppp<- cumsum(f0)*diff(x)[1]
        p0<- propTrueNullByLocalFDR(ppp)    
        f0 <- (sum(f) * f0)/sum(f0)
        cat("p0 = ",p0,"\n")
        qq[ii]<- cumsum((-f_fit$counts*log(f0/(f))) - log(p0))[ii]
    }
    qq
}
