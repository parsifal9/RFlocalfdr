#' fit.to.data.set.wrapper
#'
#' This function allows you to express your love of cats.
#' @param f_fit object returned by f.fit
#' @param imp importances
#' @param debug.flag debug flag
#' #' @keywords cats
#' @export
#' @examples
#' cat_function()
fit.to.data.set.wrapper <- function (df, imp, debug.flag = 0, plot.string = "", temp.dir = temp.dir,
                                        try.counter = 3){
    mm.df1 <- fit.to.data.set(df, imp, debug.flag = 0, plot.string = "",
                              temp.dir = temp.dir, try.counter = 1, return.all = TRUE)
    mm.df2 <- fit.to.data.set(df, imp, debug.flag = 0, plot.string = "",
                              temp.dir = temp.dir, try.counter = 2, return.all = TRUE)
    mm.df3 <- fit.to.data.set(df, imp, debug.flag = 0, plot.string = "",
                              temp.dir = temp.dir, try.counter = 3, return.all = TRUE)
    if ( class(mm.df3) == "try-error"){
        if ( class(mm.df1) == "try-error"){
            if ( class(mm.df2) == "try-error"){
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
