#' run.it.importances
#'
#' This function allows you to express your love of cats.
#' #' @param imp "reduction in impurity" importances fraom a random forest model
#' #' @param debug.flag debug.flag  either 0 (no debugging information), 1 or 2 
#' #' @keywords cats
#' @export
#' @examples
#' cat_function()
run.it.importances<-function(imp1,debug.flag=0, temp.dir=NULL,try.counter=3){

    if (debug.flag > 0){
        if ( length(temp.dir)==0){
            temp.dir <-  tempdir()
        }
        fileConn<-file(paste(temp.dir,"/output.txt",sep=""), open = "wt")
        writeLines(c("Hello","World"), fileConn)
        }
    
    imp1 <- imp1 - min(imp1) + .Machine$double.eps

    if (debug.flag > 0 ){
         writeLines(paste(length(imp1), "length(imp1)"), fileConn)
         writeLines(paste(unlist(range(imp1)), "range(imp1)"), fileConn)
    }

    f_fit<- f.fit(imp1,debug.flag=debug.flag,temp.dir=temp.dir) #makes the plot histogram_of_variable_importances.png
    y<-f_fit$zh$density
    x<-f_fit$x
    
    if (debug.flag >0){
        png(paste(temp.dir,"/density_importances.png",sep=""))
        plot(density(imp1),main="histogramand fitted spline")
        lines(x,y,col="red")
        dev.off()
    }

    df<-data.frame(x,y)
    initial.estimates <- fit.to.data.set.wraper(df,imp1,debug.flag=debug.flag,plot.string="initial",temp.dir=temp.dir,try.counter=try.counter)
    initial.estimates <-initial.estimates$Estimate

    if (debug.flag >0){
        writeLines(paste("initial estimates", initial.estimates[1],  initial.estimates[2], initial.estimates[3]), fileConn)
        writeLines("we calcualte the fdr using the initial estimates", fileConn)

        #set debug.flag=0 for local.fdr (turns off the plot  local_fdr_initial.png)
        aa<-local.fdr(f_fit,df$x,FUN=my.dsn, xi=initial.estimates[1], omega=initial.estimates[2],
                      lambda= initial.estimates[3], debug.flag=0,plot.string="initial",temp.dir=temp.dir)

        #and do our own plot with the significance level x[as.numeric(names(ww))] indicated
        png(paste(temp.dir,"/fdr_using_initial_estiamtes.png",sep=""))
        plot(x,aa,main="fdr using initial estiamtes")
        abline(h=0.2)
        ww<-which.min(abs(aa[50:119]-0.2))  #this 50 may need to be tidied up
        sum(imp1> x[as.numeric(names(ww))])  
        tt<- sum(imp1> x[as.numeric(names(ww))])
        abline(v=x[as.numeric(names(ww))])
        dev.off()

        writeLines(paste("sum(imp1> x[as.numeric(names(ww))])", tt), fileConn)
    }
    ####################################################################################################
    #take C = qsn(0.95)
    C_0.95<-sn::qsn(0.95,  xi=initial.estimates[1], omega=initial.estimates[2],  alpha= initial.estimates[3])
    if (debug.flag==1){
        writeLines(paste("calculating C_0.95", C_0.95), fileConn)
    }

    if (debug.flag >0){cat("calculating cc","\n")}
    qq <- try(
        determine.C(f_fit,df,initial.estimates,starting_value = 2,start_at=37)  ,silent = TRUE)
    if (debug.flag > 0){
        writeLines(paste(class(qq), "class(determine.C)"), fileConn)
    }
    cc <- final.estimates_cc <- aa_cc <- NA 
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


    df2<-data.frame(x[x< C_0.95],y[x< C_0.95])
    names(df2)<-c("x","y")
    final.estimates_C_0.95 <- fit.to.data.set.wrapper( df2,imp1,debug.flag=debug.flag,plot.string="final",temp.dir=temp.dir)
    #    mm1.df2.estimates <- fit.to.data.set( df2,imp1,debug.flag=debug.flag,plot.string="C",temp.dir=temp.dir)

    #should we use the cc option and C as a fallback? User settable option?
    
    if (!is.na(cc)){
        df3<-data.frame(x[x< cc],y[x< cc])
        names(df3)<-c("x","y")
        final.estimates_cc <- fit.to.data.set.wrapper( df3,imp1,debug.flag=debug.flag,plot.string="cc",temp.dir=temp.dir)
    }

    
    
    if (debug.flag > 1 ){
        #compare C and cc and the resulting fits
        png(paste(temp.dir,"/compare_C_and_cc_and_the_resulting_fits.png",sep=""))
        hist(imp1,col=6,lwd=2,breaks=100,main="compare C and cc and the resulting fits")
        abline(v=C_0.95,col="red",lwd=2)

        if (!is.na(cc)){
            abline(v= cc,col="purple",lwd=2)
        }
        legend("topright",c("C","cc"),col=c("red","purple"),lty=1)
        
        dev.off()
        
        #compare the plots
        png(paste(temp.dir,"/compare_C_and_cc_and_the_resulting_fits_2.png",sep=""))
        hist(imp1,breaks=200,freq=FALSE)
        lines(x,y,type="l",col="grey90",lwd=2,xlim=c(0,12))
        lines(df2$x,df2$y,col="green",lwd=2)
        abline(v=C_0.95,col="green")
        curve(sn::dsn(x,  xi=final.estimates_C_0.95$Estimate[1], omega=final.estimates_C_0.95$Estimate[2],
                      alpha= final.estimates_C_0.95$Estimate[3]),  add=TRUE,col="green",lwd=3)

        if (!is.na(cc)){
            lines(df3$x,df3$y,col="blue",lwd=2)
            abline(v=cc,col="blue")
            curve(sn::dsn(x,  xi=final.estimates_cc$Estimate[1], omega=final.estimates_cc$Estimate[2],
                          alpha= final.estimates_cc$Estimate[3]),  add=TRUE,col="blue",lwd=3)
        }
        legend("topright",c("C","cc"),col=c("green","blue"),lty=1)
        dev.off()
    }
    

    
    
    if (debug.flag >0){
        png(paste(temp.dir,"/fit.to.data.set_df2.png",sep=""))
        plot(x,y,type="l",col="grey90",lwd=2)
        lines(df2$x,df2$y,col="green",lwd=2)
        lines(x,my.dsn(x,xi=final.estimates_C_0.95$Estimate[1], omega=final.estimates_C_0.95$Estimate[2],  lambda= final.estimates_C_0.95$Estimate[3]))
        dev.off()
    }
    
    ###############################################################################################
    # determine p0
    ppp<-sn::psn(imp1, xi=final.estimates_C_0.95$Estimate[1], omega=final.estimates_C_0.95$Estimate[2], alpha= final.estimates_C_0.95$Estimate[3])
    p0 <- propTrueNullByLocalFDR(ppp)  
    if (debug.flag==1){
        writeLines(paste(p0, "p0"), fileConn)
    }

    ###############################################################################################
    
    
    aa_C_0.95<-local.fdr(f_fit,df$x,FUN=my.dsn, xi=final.estimates_C_0.95$Estimate[1], omega=final.estimates_C_0.95$Estimate[2],
                  lambda= final.estimates_C_0.95$Estimate[3],p0 = p0, debug.flag=debug.flag,plot.string="final",temp.dir=temp.dir)
    try(aa_cc<-local.fdr(f_fit,df$x,FUN=my.dsn, xi=final.estimates_cc$Estimate[1], omega=final.estimates_cc$Estimate[2],
                  lambda= final.estimates_cc$Estimate[3],p0 = p0, debug.flag=debug.flag,plot.string="final",temp.dir=temp.dir),silent=TRUE)
    ## mean.aa<- which.min(abs(aa-mean(imp1)))
    ## ww<-which.min(abs(aa[mean.aa:119]-0.2)) 
    ## num.sig.genes <-sum(imp1> x[as.numeric(names(ww))])  

    ## if (debug.flag >0){
    ##     png(paste(temp.dir,"/temp1.png",sep=""))
    ##     plot(x,aa_C_0.95)
    ##     abline(h=0.2)
    ##     abline(v=x[as.numeric(names(ww))])
    ##     dev.off()
    ##     cat(sum(imp1> x[as.numeric(names(ww))]),"sum(imp1> x[as.numeric(names(ww))])","\n")
    ## }
    ## if (debug.flag >0){
    ##     png(paste(temp.dir,"/temp2.png",sep=""))
    ##     hist(imp1, breaks = 200,freq=FALSE)
    ##     abline(v=sn::qsn(0.95,  xi=final.estimates_C_0.95$Estimate[1], omega=final.estimates_C_0.95$Estimate[2], alpha= final.estimates_C_0.95$Estimate[3]))
    ##     abline(v=x[as.numeric(names(ww))])
    ##     dev.off()
    ##     cat(names(imp1)[imp1> x[as.numeric(names(ww))]],"\n\n\n\n")
    ## }

    if (debug.flag >0){
        png(paste(temp.dir,"/run.it.importances_final_plot.png",sep=""))
        par(mar = c(4, 4, 4, 6)) # Set the margin on all sides to 2
        aa<- hist(imp1, col = 6, lwd = 2, breaks = 100, main = "", 
                  freq = FALSE, xlab = "importances", ylab = "density",axes=FALSE)
        curve(my.dsn(x,  xi=final.estimates_C_0.95$Estimate[1], omega=final.estimates_C_0.95$Estimate[2], lambda= final.estimates_C_0.95$Estimate[3]), 
              add = TRUE, col = "red",   lwd = 2)
        abline(v = sn::qsn(0.95, xi = final.estimates_C_0.95$Estimate[1], omega = final.estimates_C_0.95$Estimate[2], 
                            alpha = final.estimates_C_0.95$Estimate[3]), col = "red", lwd = 2)

        curve(my.dsn(x,  xi=final.estimates_cc$Estimate[1], omega=final.estimates_cc$Estimate[2], lambda= final.estimates_cc$Estimate[3]), 
              add = TRUE, col = "lawngreen"   ,   lwd = 2)
        abline(v = sn::qsn(0.95, xi = final.estimates_cc$Estimate[1], omega = final.estimates_cc$Estimate[2], 
                            alpha = final.estimates_cc$Estimate[3]), col = "lawngreen" , lwd = 2)
        
        abline(v = x[as.numeric(names(ww))], lwd = 2,  col = "orange")
        #looks like it misses the intersection
        abline(v = C_0.95, lwd = 2, col = "blue")
        abline(v = cc, lwd = 2, col = "purple")
        
        axis(2, pretty( c(0,max(aa$density)+0.5*max(aa$density)),10))
        par(new=TRUE)
        plot(c(0,max(x)),c(0,1),type="n",axes=FALSE, xlab = "", ylab = "")
        lines(x, aa_C_0.95, lwd = 3)
        abline(h = 0.2,col="orange")
        axis(4, pretty(c(0,1),10))
        mtext(side=4,line=3,  "local fdr")
        axis(1,pretty(range(1:max(x)),10))
        legend("topright", c("fitted curve", "95% quantile", 
                             "cutoff", "C", "cc","fdr"), col = c("red", "red", "orange", 
                                                                 "blue", "purple","black"), lty = 1, lwd = 2)
        box() #
        dev.off()
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

    temp<-list(aa_C_0.95,aa_cc,df$x,final.estimates_C_0.95,final.estimates_cc,temp.dir,C_0.95,cc,p0)
    names(temp)<-c("fdr_0.95","fdr_cc","x","estimates_C_0.95,","estimates_cc","temp.dir","C_0.95","cc","p0")
    temp
}










