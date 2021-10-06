#' run.it.importances
#'
#' This function allows you to express your love of cats.
#' #' @param imp importances
#' #' @param debyg debug
#' #' @keywords cats
#' @export
#' @examples
#' cat_function()
run.it.importances<-function(imp1,debug=0){
    imp1 <- imp1 - min(imp1) + .Machine$double.eps
    if (debug==1){
        cat(length(imp1), "length(imp1)","\n")
        cat(range(imp1), "range(imp1)","\n")
        plot(density(imp1))
    }

    f_fit<- f.fit(imp1)
    y<-f_fit$zh$density
    x<-f_fit$x
    if (debug==1){
        plot(x,y,type="l")
    }

    df<-data.frame(x,y)
    initial.estimates <- fit.to.data.set(df,imp1)
    initial.estimates <-initial.estimates$Estimate

    if (debug==1){
        cat("initial estimates", initial.estimates[1],  initial.estimates[2], initial.estimates[3],"\n")
    }

    if (debug==1){
        aa<-local.fdr(f_fit,df$x,FUN=my.dsn, xi=initial.estimates[1], omega=initial.estimates[2],
                      lambda= initial.estimates[3])
        plot(x,aa)
        abline(h=0.2)
        ww<-which.min(abs(aa[50:119]-0.2))
        sum(imp1> x[as.numeric(names(ww))])  
        cat( sum(imp1> x[as.numeric(names(ww))]), " sum(imp1> x[as.numeric(names(ww))])","\n")
        abline(v=x[as.numeric(names(ww))])
    }
    ####################################################################################################
    if (debug==1){cat("calculating C","\n")}
    #take C = qsn(0.95)
    C<-sn::qsn(0.95,  xi=initial.estimates[1], omega=initial.estimates[2],  alpha= initial.estimates[3])
    if (debug==1){cat(C, "C","\n")}

    if (debug==1){cat("calculating cc","\n")}
    qq<-try(
        determine.C(f_fit,df,initial.estimates,starting_value = 2,start_at=37)  ,silent = TRUE)
    if (debug==1){cat(class(qq), "class(determine.C) ","\n")}
    if (class(qq)=="numeric"){
        cc<-x[which.min(qq)] 
        cat(cc, "cc","\n")
        plot(x,qq)
        abline(v=cc)
   }

    if (FALSE){
        hist(imp1,col=6,lwd=2,breaks=100,main="histogram of importances")
        abline(v=C,col="red",lwd=2)
        abline(v= cc,col="purple",lwd=2)
        
        
        df2<-data.frame(x[x< C],y[x< C])
        names(df2)<-c("x","y")
        
        mm1.df2 = minpack.lm::nlsLM(y ~ my.dsn(x,xi=xi, omega=omega, lambda=lambda), 
                                    start = list(xi=mean(x), omega=2, lambda= 1 ),
                                    data = df2,control=minpack.lm::nls.lm.control(maxiter=400))
        
        df3<-data.frame(x[x< cc],y[x< cc])
        names(df2)<-c("x","y")
        
        mm1.df3 = minpack.lm::nlsLM(y ~ my.dsn(x,xi=xi, omega=omega, lambda=lambda), 
                                    start = list(xi=mean(x), omega=2, lambda= 1 ),
                                    data = df3,control=minpack.lm::nls.lm.control(maxiter=400))
        
        
        #compare the plots 
        plot(x,y,type="l",col="grey90",lwd=2,xlim=c(0,12))
        lines(df2$x,df2$y,col="green",lwd=2)
        lines(x,predict(mm1.df2,newdata=df),col="green",lwd=3)  
        
        lines(df3$x,df3$y,col="blue",lwd=2)
        lines(x,predict(mm1.df3,newdata=df),col="blue",lwd=3)  
    }
    
    final.estimates <- fit.to.data.set( df2<-data.frame(x[x< C],y[x< C]),imp1)$Estimate

    if (debug==1){
     plot(x,y,type="l",col="grey90",lwd=2,xlim=c(0,12))
     lines(df2$x,df2$y,col="green",lwd=2)
     lines(x,my.dsn(x,xi=final.estimates[1], omega=final.estimates[2],
                               lambda= final.estimates[3]))
    }
    
    ###############################################################################################
    # determine p0
    ppp<-sn::psn(imp1, xi=final.estimates[1], omega=final.estimates[2], alpha= final.estimates[3])
    p0<- propTrueNullByLocalFDR(ppp)  
    if (debug==1){
        cat(p0, "p0","\n")
    }

    ###############################################################################################
    
    
    aa<-local.fdr(f_fit,df$x,FUN=my.dsn, xi=final.estimates[1], omega=final.estimates[2],
                  lambda= final.estimates[3],p0 = p0 )
    ## mean.aa<- which.min(abs(aa-mean(imp1)))
    ## ww<-which.min(abs(aa[mean.aa:119]-0.2)) 
    ## num.sig.genes <-sum(imp1> x[as.numeric(names(ww))])  

    ## if (debug==1){
    ## plot(x,aa)
    ## abline(h=0.2)
    ## abline(v=x[as.numeric(names(ww))])
    ## cat(sum(imp1> x[as.numeric(names(ww))]),"sum(imp1> x[as.numeric(names(ww))])","\n")
    ## hist(imp1, breaks = 200,freq=FALSE)
    ## abline(v=sn::qsn(0.95,  xi=final.estimates[1], omega=final.estimates[2], alpha= final.estimates[3]))
    ## abline(v=x[as.numeric(names(ww))])
    ## cat(names(imp1)[imp1> x[as.numeric(names(ww))]],"\n\n\n\n")
    ## }
    
    ## a1<- match(names(imp1)[imp1> x[as.numeric(names(ww))]],names(imp1))
    ## ppp<- 1-sn::psn(imp1[a1], xi=final.estimates[1], omega=final.estimates[2], alpha= final.estimates[3])
    ## names(ppp) <-names(imp1)[imp1> x[as.numeric(names(ww))]]
    ## ppp
    temp<-list(aa,df$x,final.estimates)
    names(temp)<-c("fdr","x","estimates")
    temp
}










