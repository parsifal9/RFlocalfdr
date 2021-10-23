#' run.it.importances
#'
#' This function allows you to express your love of cats.
#' #' @param imp importances
#' #' @param debug.flag debug
#' #' @keywords cats
#' @export
#' @examples
#' cat_function()
run.it.importances<-function(imp1,debug.flag=0){
    temp.dir <- NULL

    if (debug.flag > 0){
        temp.dir <-  tempdir()
        fileConn<-file(paste(temp.dir,"/output.txt",sep=""), open = "wt")
        writeLines(c("Hello","World"), fileConn)
        }
    
    imp1 <- imp1 - min(imp1) + .Machine$double.eps

    if (debug.flag > 0 ){
         writeLines(paste(length(imp1), "length(imp1)"), fileConn)
         writeLines(paste(unlist(range(imp1)), "range(imp1)"), fileConn)
         png(paste(temp.dir,"/density_importances.png",sep=""))
         plot(density(imp1))
         dev.off()
    }

    f_fit<- f.fit(imp1,debug.flag=debug.flag,temp.dir=temp.dir)
    y<-f_fit$zh$density
    x<-f_fit$x
    
    if (debug.flag >0){
        png(paste(temp.dir,"/density_importances2.png",sep=""))
        plot(x,y,type="l")
        dev.off()
    }

    df<-data.frame(x,y)
    initial.estimates <- fit.to.data.set(df,imp1,debug.flag=debug.flag,plot.string="initial",temp.dir=temp.dir)
    initial.estimates <-initial.estimates$Estimate

    if (debug.flag >0){
        writeLines(paste("initial estimates", initial.estimates[1],  initial.estimates[2], initial.estimates[3]), fileConn)

        writeLines("we calcualte the fdr using the initial estimates", fileConn)
        aa<-local.fdr(f_fit,df$x,FUN=my.dsn, xi=initial.estimates[1], omega=initial.estimates[2],
                      lambda= initial.estimates[3])  #return a plot  -- do we want this? 
        plot(x,aa,main="fdr using initial estiamtes")
        abline(h=0.2)
        ww<-which.min(abs(aa[50:119]-0.2))  #this 50 may need to be tidied up
        sum(imp1> x[as.numeric(names(ww))])  
        tt<- sum(imp1> x[as.numeric(names(ww))])
        abline(v=x[as.numeric(names(ww))])
        writeLines(paste("sum(imp1> x[as.numeric(names(ww))])", tt), fileConn)
    }
    ####################################################################################################
    #take C = qsn(0.95)
    C<-sn::qsn(0.95,  xi=initial.estimates[1], omega=initial.estimates[2],  alpha= initial.estimates[3])
    if (debug.flag==1){
        writeLines(paste("calculating C", C), fileConn)
    }

    if (debug.flag >0){cat("calculating cc","\n")}
    qq <- try(
        determine.C(f_fit,df,initial.estimates,starting_value = 2,start_at=37)  ,silent = TRUE)
    if (debug.flag > 0){
        writeLines(paste(class(qq), "class(determine.C)"), fileConn)
    }
    cc <- NA
    if (class(qq) == "numeric"){
        cc<-x[which.min(qq)]
        if (debug.flag > 0 ){
            writeLines(paste("cc= ",cc), fileConn)
            png(paste(temp.dir,"/determine_cc.png",sep=""))
            plot(x,qq,main="determine cc")
            abline(v=cc)
            dev.off()
        }

    }

    if (debug.flag > 1 ){
        #compare C and cc and the resulting fits
        png(paste(temp.dir,"/compare_C_and_cc_and_the_resulting_fits.png",sep=""))
        hist(imp1,col=6,lwd=2,breaks=100,main="compare C and cc and the resulting fits")
        abline(v=C,col="red",lwd=2)

        if (!is.na(cc)){
            abline(v= cc,col="purple",lwd=2)
        }

        dev.off()

        
        df2<-data.frame(x[x< C],y[x< C])
        names(df2)<-c("x","y")
        mm1.df2.estimates <- fit.to.data.set( df2,imp1,debug.flag=debug.flag,plot.string="C",temp.dir=temp.dir)

        if (!is.na(cc)){
            df3<-data.frame(x[x< cc],y[x< cc])
            names(df3)<-c("x","y")
            mm1.df3.estimates <- fit.to.data.set( df3,imp1,debug.flag=debug.flag,plot.string="cc",temp.dir=temp.dir)
         }
        
        #compare the plots
        hist(imp,breaks=200,freq=FALSE)
        lines(x,y,type="l",col="grey90",lwd=2,xlim=c(0,12))
        lines(df2$x,df2$y,col="green",lwd=2)
        abline(v=C,col="green")
        curve(sn::dsn(x,  xi=mm1.df2.estimates$Estimate[1], omega=mm1.df2.estimates$Estimate[2],
                      alpha= mm1.df2.estimates$Estimate[3]),  add=TRUE,col="green",lwd=3)

        if (!is.na(cc)){
            lines(df3$x,df3$y,col="blue",lwd=2)
            abline(v=cc,col="blue")
            curve(sn::dsn(x,  xi=mm1.df3.estimates$Estimate[1], omega=mm1.df3.estimates$Estimate[2],
                          alpha= mm1.df3.estimates$Estimate[3]),  add=TRUE,col="blue",lwd=3)
        }
    }
    
    final.estimates <- fit.to.data.set( df2<-data.frame(x[x< C],y[x< C]),imp1,debug.flag=debug.flag,plot.string="final",temp.dir=temp.dir)$Estimate
    #should we use the cc option and C as a fallback? Usersettable option?

    
    if (debug.flag >0){
        png(paste(temp.dir,"/fit.to.data.set_df2.png",sep=""))
        plot(x,y,type="l",col="grey90",lwd=2)
        lines(df2$x,df2$y,col="green",lwd=2)
        lines(x,my.dsn(x,xi=final.estimates[1], omega=final.estimates[2],
                       lambda= final.estimates[3]))
        dev.off()
    }
    
    ###############################################################################################
    # determine p0
    ppp<-sn::psn(imp1, xi=final.estimates[1], omega=final.estimates[2], alpha= final.estimates[3])
    p0 <- propTrueNullByLocalFDR(ppp)  
    if (debug.flag==1){
        writeLines(paste(p0, "p0"), fileConn)
    }

    ###############################################################################################
    
    
    aa<-local.fdr(f_fit,df$x,FUN=my.dsn, xi=final.estimates[1], omega=final.estimates[2],
                  lambda= final.estimates[3],p0 = p0 )
    ## mean.aa<- which.min(abs(aa-mean(imp1)))
    ## ww<-which.min(abs(aa[mean.aa:119]-0.2)) 
    ## num.sig.genes <-sum(imp1> x[as.numeric(names(ww))])  

    if (debug.flag >0){
        plot(x,aa)
        abline(h=0.2)
        abline(v=x[as.numeric(names(ww))])
        cat(sum(imp1> x[as.numeric(names(ww))]),"sum(imp1> x[as.numeric(names(ww))])","\n")
    }
    if (debug.flag >0){
        hist(imp1, breaks = 200,freq=FALSE)
        abline(v=sn::qsn(0.95,  xi=final.estimates[1], omega=final.estimates[2], alpha= final.estimates[3]))
        abline(v=x[as.numeric(names(ww))])
        cat(names(imp1)[imp1> x[as.numeric(names(ww))]],"\n\n\n\n")
    }


    if (debug.flag >1){
        cat(names(imp1)[imp1> x[as.numeric(names(ww))]],"\n\n\n\n")
    }

    ## a1<- match(names(imp1)[imp1> x[as.numeric(names(ww))]],names(imp1))
    ## ppp<- 1-sn::psn(imp1[a1], xi=final.estimates[1], omega=final.estimates[2], alpha= final.estimates[3])
    ## names(ppp) <-names(imp1)[imp1> x[as.numeric(names(ww))]]
    ## ppp

    if (debug.flag > 0){
        close(fileConn)
    }

    temp<-list(aa,df$x,final.estimates,temp.dir)
    names(temp)<-c("fdr","x","estimates","temp.dir")
    temp
}










