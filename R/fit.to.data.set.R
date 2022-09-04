#' fit.to.data.set
#'
#' This function allows you to express your love of cats.
#' @param f_fit object returned by f.fit
#' @param imp importances
#' @param debug.flag debug flag
#' #' @keywords cats
#' @export
#' @examples
#' cat_function()
fit.to.data.set<-function(df,imp,debug.flag=0,plot.string="",temp.dir=NULL,try.counter=3,return.all=FALSE){

    mm1.df <- NA
    class(mm1.df) <- "try-error"
    x<-df$x
    y<-df$y

    if (class(mm1.df)=="try-error" & try.counter == 0 ){
        mm1.df <-try(
            minpack.lm::nlsLM(y ~ my.dsn(x,xi=xi, omega=omega, lambda=lambda), 
                              start = list(xi=  1   , omega=2, lambda= 1   ),
                              data = df,control=minpack.lm::nls.lm.control(maxiter=400))
           ,silent = TRUE)
        try.counter <- 1
        
        
        if  (debug.flag >1 ){
            cat(class(mm1.df), "class(mm1.df) -- try 1","\n")
        }
    }

    if (class(mm1.df)=="try-error" & (try.counter < 2)){
        mm1.df <-try(
            minpack.lm::nlsLM(y ~ my.dsn(x,xi=xi, omega=omega, lambda=lambda), 
                              start = list(xi=  mean(x)   , omega=2, lambda= 1   ),
                              data = df,control=minpack.lm::nls.lm.control(maxiter=400))
           ,silent = TRUE)
        try.counter <- 2
        
        
        if  (debug.flag >1 ){
            cat(class(mm1.df), "class(mm1.df) -- try 2","\n")
        }
    }

    if (class(mm1.df) == "try-error" & try.counter == 3){
        vip.sn.mle <- fitdistrplus::fitdist(imp, "sn",start=list( xi = mean(imp),omega = 1, alpha= 0),
                              lower=c(-Inf,0,-Inf),fix.arg=list( tau=0),method = c("mle"))
        mm1.df <-try( minpack.lm::nlsLM(y ~ my.dsn(x,xi=xi, omega=omega, lambda=lambda), 
                                        start = list(xi=  vip.sn.mle$estimate[1], omega=vip.sn.mle$estimate[2],
                                                     lambda= vip.sn.mle$estimate[3]),
                                        data = df,control=minpack.lm::nls.lm.control(maxiter=400)) ,silent = TRUE)
        try.counter <- 3
    }
    
    if  (debug.flag >1 ){
    cat(class(mm1.df), "class(mm1.df)-- try 3","\n")
    }
    
    if (debug.flag >0 ){
        cat("try =",try.counter, "\n")
    }

    if (debug.flag >1 ){
        png(paste(temp.dir,"/fit_to_data_set_",plot.string,".png",sep=""))
        hist(imp,breaks=200,freq=FALSE)
        lines(df$x,df$y,type="l",col="green",lwd=2,xlim=c(0, max(df$x)+0.5))
      

        if (try.counter ==3){
        curve(sn::dsn(x, xi= vip.sn.mle$estimate[1], omega=vip.sn.mle$estimate[2], alpha= vip.sn.mle$estimate[3]),
              add=TRUE,col="purple",lwd=3,xlim=c(0,16))
        curve(my.dsn(x,xi=  vip.sn.mle$estimate[1], omega=vip.sn.mle$estimate[2], lambda= vip.sn.mle$estimate[3]),
              add=TRUE,col="purple",lwd=3)

        }
        
        lines(x,predict(mm1.df),col="red",lwd=3)
        curve(my.dsn(x,  xi=summary(mm1.df)$parameters[1,1], omega=summary(mm1.df)$parameters[2,1],
                     lambda= summary(mm1.df)$parameters[3,1] ),  add=TRUE,col="blue",lwd=3)
        curve(sn::dsn(x,  xi=summary(mm1.df)$parameters[1,1], omega=summary(mm1.df)$parameters[2,1],
                      alpha= summary(mm1.df)$parameters[3,1] ),  add=TRUE,col="blue",lwd=3)

        abline(v=sn::qsn(0.5,  xi=summary(mm1.df)$parameters[1,1], omega=summary(mm1.df)$parameters[2,1],
                         alpha= summary(mm1.df)$parameters[3,1] ),col="gray30")

        try(abline(v=sn::qsn(0.05,  xi=summary(mm1.df)$parameters[1,1], omega=summary(mm1.df)$parameters[2,1],
                         alpha= summary(mm1.df)$parameters[3,1] )),silent=TRUE)

        abline(v=sn::qsn(0.95,  xi=summary(mm1.df)$parameters[1,1], omega=summary(mm1.df)$parameters[2,1],
                         alpha= summary(mm1.df)$parameters[3,1] ),col="gray30")
        
        cat(sum(abs(df$y-predict(mm1.df))), "sum(abs(df$y-predict(mm1.df)))","\n")
        
        if (try.counter ==3){
            legend("topright",
                   c("spline fit","fitted f0","initial fitdist fit", "Skew-normal at fitted values","quantiles of skew normal"), col=c("green","red","purple","blue","grey30"),
                   lty=1, lwd=4)
        }else{
            legend("topright",
                   c("spline fit","fitted f0","Skew-normal at fitted values","quantiles of skew normal"), col=c("green","red","blue","grey30"),
                   lty=1, lwd=4)
        }
        dev.off()
    }

    if (return.all=TRUE){
        aa<-mm1.df
    } else{
        aa<- data.frame(summary(mm1.df)$parameters)
    }
    aa
}
