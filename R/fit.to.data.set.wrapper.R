#' fit.to.data.set.wrapper
#'
#' This function allows you to express your love of cats.
#' @param imp importances
#' @param debug.flag debug flag
#' @param df, contains x and y, midpoints and counts from a histogram of imp
#' @param plot.string, file name for a debugging plot, passed on to fit.to.data.set
#' @param temp.dir, directory for debugging output, passed on to fit.to.data.set
#' @param try.counter, passed on to fit.to.data.set
#'        try.counter=1 my.dsn xi=  1  
#'        try.counter=2  xi=  mean(x)
#'        try.counter=3 start xi, omega, lambda from the parameters retuned by fitdistrplus::fitdist
#' @export
#' @return  If the skew-normal fitting routine is succesful, then the matrix of parmaters and standard errors is returned.
#'         -- othewise a "try-error" message is returned
#' @examples
#' \dontrun{
#' data(ch22)
#' t2 <-ch22$C
#' imp<-log(ch22$imp)
#' imp<-imp[t2 > 30]
#' imp <- imp[imp != -Inf]
#' imp <- imp - min(imp) + .Machine$double.eps
#' f_fit <- f.fit(imp )
#' y <- f_fit$zh$density
#' x <- f_fit$midpoints
#' C <- quantile(imp,probs=0.75)
#' df2 <- data.frame(x[x < C], y[x < C])
#' initial.estimates <- fit.to.data.set.wrapper(df2, imp)
#' #Nonlinear regression model                                            
#' #  model: y ~ my.dsn(x, xi = xi, omega = omega, lambda = lambda)       
#' #   data: df                                                           
#' #       xi.xi  omega.omega lambda.alpha                                
#' #       1.163        1.193        1.574                                
#' # residual sum-of-squares: 0.06269                                     
#' #                                                                      
#' #Number of iterations to convergence: 23                               
#' #Achieved convergence tolerance: 1.49e-08
#' }
  


fit.to.data.set.wrapper <- function (df, imp, debug.flag = 0, plot.string = "", temp.dir = temp.dir,
                                     return.all = TRUE,
                                     try.counter = 3){

    if (any(!is.finite(imp))){
        stop("importance is not finite")
    }
    if (try.counter == 1){
    mm.df1 <- fit.to.data.set(df, imp, debug.flag = debug.flag, plot.string = plot.string,
                              temp.dir = temp.dir, try.counter = 1, return.all =  return.all)
      return( mm.df1)
    }
    if (try.counter == 2){
    mm.df2 <- fit.to.data.set(df, imp, debug.flag = debug.flag, plot.string = plot.string,
                              temp.dir = temp.dir, try.counter = 2, return.all =  return.all)
     return( mm.df2)
    }
    if (try.counter == 3){
    mm.df3 <- fit.to.data.set(df, imp, debug.flag = debug.flag, plot.string = plot.string,
                              temp.dir = temp.dir, try.counter = 3, return.all =  return.all)
    return( mm.df3)
    }

    if (try.counter == "all"){
        mm.df1 <- fit.to.data.set(df, imp, debug.flag = debug.flag, plot.string = plot.string,
                                  temp.dir = temp.dir, try.counter = 1, return.all =  return.all)
        mm.df2 <- fit.to.data.set(df, imp, debug.flag = debug.flag, plot.string = "",
                                  temp.dir = temp.dir, try.counter = 2, return.all =  return.all)
        mm.df3 <- fit.to.data.set(df, imp, debug.flag = debug.flag, plot.string = "",
                                  temp.dir = temp.dir, try.counter = 3, return.all =  return.all)
        
        if ( inherits(mm.df3,"try-error")){
            if ( inherits(mm.df1,"try-error")){
                if ( inherits(mm.df2,"try-error")){
                    return("try-error")
                } else {
                    return(  mm.df2)
                }
            }
            else{ return( mm.df1)}
        } else{
            return( mm.df3)
        }

    }
}

## fit.to.data.set.wrapper<-function(df,imp,debug.flag=0,plot.string="",temp.dir=temp.dir,try.counter=3){

##     mm.df1<-fit.to.data.set(df,imp,debug.flag=0,plot.string="",temp.dir=temp.dir,try.counter=1,return.all=TRUE)
##     mm.df2<-fit.to.data.set(df,imp,debug.flag=0,plot.string="",temp.dir=temp.dir,try.counter=2,return.all=TRUE)
##     mm.df3<-fit.to.data.set(df,imp,debug.flag=0,plot.string="",temp.dir=temp.dir,try.counter=3,return.all=TRUE)

##     mm.df3

## }
