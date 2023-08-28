#' fit.to.data.set.wrapper
#'
#' This function allows you to express your love of cats.
#' @param imp importances
#' @param debug.flag debug flag
#' @param df, contains x and y, midpoints and counts from a histogram of imp
#' @param plot.string, file name for a debugging plot
#' @param temp.dir, directory for debugging output
#' @param try.counter
#'        try.counter=1 my.dsn xi=  1  
#'        try.counter=2  xi=  mean(x)
#'        try.counter=3 start xi, omega, lambda from the parameters retuned by fitdistrplus::fitdist
#' @export
#' @examples
#' \dontrun{
#' data(ch22)
#' t2 <-ch22$C
#' imp<-log(ch22$imp)
#' temp<-imp[t2 > 30]
#' temp <- temp[temp != -Inf]
#' temp <- temp - min(temp) + .Machine$double.eps
#' f_fit <- f.fit(temp )
#' y <- f_fit$zh$density
#' x <- f_fit$midpoints
#' C <- quantile(temp,probs=0.75)
#' df2 <- data.frame(x[x < C], y[x < C])
#' initial.estimates <- fit.to.data.set.wrapper(df2, temp)
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
                                        try.counter = 3){
    mm.df1 <- fit.to.data.set(df, imp, debug.flag = 0, plot.string = "",
                              temp.dir = temp.dir, try.counter = 1, return.all = TRUE)
    mm.df2 <- fit.to.data.set(df, imp, debug.flag = 0, plot.string = "",
                              temp.dir = temp.dir, try.counter = 2, return.all = TRUE)
    mm.df3 <- fit.to.data.set(df, imp, debug.flag = 0, plot.string = "",
                              temp.dir = temp.dir, try.counter = 3, return.all = TRUE)
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

## fit.to.data.set.wrapper<-function(df,imp,debug.flag=0,plot.string="",temp.dir=temp.dir,try.counter=3){

##     mm.df1<-fit.to.data.set(df,imp,debug.flag=0,plot.string="",temp.dir=temp.dir,try.counter=1,return.all=TRUE)
##     mm.df2<-fit.to.data.set(df,imp,debug.flag=0,plot.string="",temp.dir=temp.dir,try.counter=2,return.all=TRUE)
##     mm.df3<-fit.to.data.set(df,imp,debug.flag=0,plot.string="",temp.dir=temp.dir,try.counter=3,return.all=TRUE)

##     mm.df3

## }
