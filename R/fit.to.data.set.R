#' fit.to.data.set
#'
#' This function allows you to express your love of cats.
#' @param f_fit object returned by f.fit
#' @param imp importances
#' #' @keywords cats
#' @export
#' @examples
#' cat_function()
fit.to.data.set<-function(df,imp){
    #just a test fit to the whole data set
    #what about imp1 ?????????????
    x<-df$x
    y<-df$y
    mm1.df <-try(
        minpack.lm::nlsLM(y ~ my.dsn(x,xi=xi, omega=omega, lambda=lambda), 
                          start = list(xi=  1   , omega=2, lambda= 1   ),
                          data = df,control=minpack.lm::nls.lm.control(maxiter=400))
       ,silent = TRUE)
    
    cat(class(mm1.df), "class(mm1.df) -- try 1","\n")
    
    if (class(mm1.df)=="try-error"){
        mm1.df <-try(
            minpack.lm::nlsLM(y ~ my.dsn(x,xi=xi, omega=omega, lambda=lambda), 
                              start = list(xi=  mean(x)   , omega=2, lambda= 1   ),
                              data = df,control=minpack.lm::nls.lm.control(maxiter=400))
           ,silent = TRUE)
    }
    
    cat(class(mm1.df), "class(mm1.df) -- try 2","\n")
    
    if (class(mm1.df)=="try-error"){
        vip.sn.mle <- fitdistrplus::fitdist(imp, "sn",start=list( xi = mean(imp),omega = 1, alpha= 0),
                              lower=c(-Inf,0,-Inf),fix.arg=list( tau=0),method = c("mle"))
        plot(df$x,df$y,type="l",col="green",lwd=2,xlim=c(0, max(df$x)+0.5))
        curve(sn::dsn(x, xi= vip.sn.mle$estimate[1], omega=vip.sn.mle$estimate[2], alpha= vip.sn.mle$estimate[3]),
              add=TRUE,col="blue",lwd=3,xlim=c(0,16))
        curve(my.dsn(x,xi=  vip.sn.mle$estimate[1], omega=vip.sn.mle$estimate[2], lambda= vip.sn.mle$estimate[3]),
              add=TRUE,col="purple",lwd=3)
        mm1.df <-try( minpack.lm::nlsLM(y ~ my.dsn(x,xi=xi, omega=omega, lambda=lambda), 
                                        start = list(xi=  vip.sn.mle$estimate[1], omega=vip.sn.mle$estimate[2],
                                                     lambda= vip.sn.mle$estimate[3]),
                                        data = df,control=minpack.lm::nls.lm.control(maxiter=400)) ,silent = TRUE)
    }
    
    if (FALSE){
        ## mm1.df <-  minpack.lm::nlsLM(y ~ my.dsn(x,xi=xi, omega=omega, lambda=lambda), 
        ##                              start = list(xi=  10, omega=1, lambda= 1),
        ##                              data = df,control=minpack.lm::nls.lm.control(maxiter=400)) #
        
        plot(df$x,df$y,type="l",col="green",lwd=2,xlim=c(0, max(df$x)+0.5))
        curve(my.dsn(x,xi=  vip.sn.mle$estimate[1], omega=vip.sn.mle$estimate[2], lambda= vip.sn.mle$estimate[3]),
              add=TRUE,col="purple",lwd=3)
        curve(sn::dsn(x, xi= vip.sn.mle$estimate[1], omega=vip.sn.mle$estimate[2], alpha= vip.sn.mle$estimate[3]),
              add=TRUE,col="blue",lwd=3)
        lines(x,predict(mm1.df),col="red",lwd=3)
        system("sleep 5")
    }
    
    cat(class(mm1.df), "class(mm1.df)-- try 3","\n")
    if (FALSE){
        plot(df$x,df$y,type="l",col="green",lwd=2,xlim=c(0, max(df$x)+0.5))
        lines(x,predict(mm1.df),col="red",lwd=3)  
        curve(my.dsn(x,  xi=summary(mm1.df)$parameters[1,1], omega=summary(mm1.df)$parameters[2,1],
                     lambda= summary(mm1.df)$parameters[3,1] ),  add=TRUE,col="blue",lwd=3)
        curve(sn::dsn(x,  xi=summary(mm1.df)$parameters[1,1], omega=summary(mm1.df)$parameters[2,1],
                      alpha= summary(mm1.df)$parameters[3,1] ),  add=TRUE,col="blue",lwd=3)
        abline(v=sn::qsn(0.5,  xi=summary(mm1.df)$parameters[1,1], omega=summary(mm1.df)$parameters[2,1],
                         alpha= summary(mm1.df)$parameters[3,1] ))
        abline(v=sn::qsn(0.05,  xi=summary(mm1.df)$parameters[1,1], omega=summary(mm1.df)$parameters[2,1],
                         alpha= summary(mm1.df)$parameters[3,1] ))
        abline(v=sn::qsn(0.95,  xi=summary(mm1.df)$parameters[1,1], omega=summary(mm1.df)$parameters[2,1],
                         alpha= summary(mm1.df)$parameters[3,1] ))
        sum(abs(df$y-predict(mm1.df))) # 
        cat(sum(abs(df$y-predict(mm1.df))), "sum(abs(df$y-predict(mm1.df)))","\n")
        }
    data.frame(summary(mm1.df)$parameters)
}
