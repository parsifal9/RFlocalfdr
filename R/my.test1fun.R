#' my.test1fun
#'
#' tests the compliance of skew-normal distribution functions with the expectations
#' of the package fitdistrplus
#'
#' RFlocalfdr also uses wrappers around the funcions dsn, qsn and psn from the package "sn"
#' https://cran.r-project.org/web/packages/sn/.
#' This is due to the fact that 
#' fitdistrplus::fitdist(imp, "sn", start = list(xi = mean(imp)...
#' returns warnings such as
#' The dsn function should return a zero-length vector when input has length zero and not raise an error
#' The psn function should have its first argument named: q as in base R
#' These wrappers ensure conformity with the expectations of fitdistrplus::fitdist
#' 
#' @param fn - name of the function to be tested "dsn", "psn" or "qsn"
#' @param start.arg - the starting arguments for fitting the density
#' @param fix.arg - fixed arguments
#' @param dpqr -- are we testing the "d", "p" or "q" function? not needed as it can be inferred from the argument "fn"
#' @keywords skew normal
#' @export my.test1fun
#' @return fits the Density function for the skew-normal (SN) distribution.
#' @examples
#' library(sn)
#' curve(sn::dsn(x,xi=0, omega=1, alpha=1, tau=0),xlim=c(-10,10),col="blue")
#' curve(sn::dsn(x,xi=0, omega=1, alpha=0.1, tau=0),xlim=c(-10,10),col="blue",add=TRUE)
#' curve(sn::dsn(x,xi=1, omega=2, alpha=2, tau=0),xlim=c(-10,20),col="blue",add=TRUE)
#' curve(sn::dsn(x,xi=3, omega=4, alpha=4, tau=0),xlim=c(-10,20),col="blue",add=TRUE)
#'
#' curve(my.dsn(x),xlim=c(-10,10),col="red",add=TRUE)
#' curve(my.dsn(x,lambda=0.1),xlim=c(-10,10),col="red",add=TRUE)
#' curve(my.dsn(x,xi=1, omega=2, lambda=2),xlim=c(-10,20),col="red",add=TRUE)
#' curve(my.dsn(x,xi=3, omega=4, lambda=4),xlim=c(-10,20),col="red",add=TRUE)
#'
#' #dsn, qsn and psn are wrappers around the provided functions provided by sn. This is done to
#' # overcome some checking done by fitdistrplus
#' \donttest{
#' library(sn)
#' getAnywhere("dsn")
#' RFlocalfdr::my.test1fun("sn::dsn", list(xi = -Inf, omega =1, alpha=0 ), fix.arg = list(tau = 0))
#' RFlocalfdr::my.test1fun("sn::psn", list(xi = -Inf, omega =1, alpha=0 ), fix.arg = list(tau = 0))
#' RFlocalfdr::my.test1fun("sn::qsn", list(xi = -Inf, omega =1, alpha=0 ), fix.arg = list(tau = 0))
#' #all return FALSE
#'
#' detach("package:sn", unload=TRUE)
#' getAnywhere("dsn")
#' RFlocalfdr::my.test1fun("dsn", list(xi = -Inf, omega =1, alpha=0 ), fix.arg = list(tau = 0))#TRUE
#' RFlocalfdr::my.test1fun("psn", list(xi = -Inf, omega =1, alpha=0 ), fix.arg = list(tau = 0))#TRUE
#' RFlocalfdr::my.test1fun("qsn", list(xi = -Inf, omega =1, alpha=0 ), fix.arg = list(tau = 0))#TRUE
#' }



my.test1fun<- function (fn, start.arg, fix.arg, dpqr) 
{
    res <- data.frame(ok = FALSE, txt = "")
    stopifnot(is.list(start.arg))
    if (!is.null(fix.arg)) 
        stopifnot(is.list(fix.arg))
    if (!exists(fn, mode = "function")) {
        res$txt <- paste("The", fn, "function must be defined")
        return(res)
    }
    if (missing(dpqr)) 
        dpqr <- substr(fn, 1, 1)
    firstarg_theo <- switch(dpqr, d = "x", p = "q", q = "p",   r = "n")
    firstarg_found <- names(formals(fn))[1]
    if (firstarg_found != firstarg_theo) {
        t0 <- paste("The", fn, "function should have its first argument named:", 
            firstarg_theo)
        res$txt <- paste(t0, "as in base R")
        return(res)
    }
    res0 <- try(do.call(fn, c(list(numeric(0)), start.arg, fix.arg)),    silent = TRUE)
    t0 <- paste("The", fn, "function should return a zero-length vector when input has length zero and not raise an error")
    t1 <- paste("The", fn, "function should return a zero-length vector when input has length zero")
    if (inherits(res0, "try-error")) {
        res$txt <- t0
        return(res)
    }
    if (length(res0) != 0) {
        res$txt <- t1
        return(res)
    }
    x <- c(0, 1, Inf, NaN, -1)
    res1 <- try(do.call(fn, c(list(x), start.arg, fix.arg)), silent = TRUE)
    t2 <- paste("The", fn, "function should return a vector of with NaN values when input has inconsistent values and not raise an error")
    if (inherits(res1, "try-error")) {
        res$txt <- t2
        return(res)
    }
    x <- c(0, 1, NA)
    res2 <- try(do.call(fn, c(list(x), start.arg, fix.arg)),   silent = TRUE)
    t4 <- paste("The", fn, "function should return a vector of with NA values when input has missing values and not raise an error")
    t5 <- paste("The", fn, "function should return a vector of with NA values when input has missing values and not remove missing values")
    if (inherits(res2, "try-error")) {
        res$txt <- t4
        return(res)
    }
    if (length(res2) != length(x)) {
        res$txt <- t5
        return(res)
    }
    x <- 0:1
    start.arg <- lapply(start.arg, function(x) -x)
    res3 <- try(do.call(fn, c(list(x), start.arg, fix.arg)),   silent = TRUE)
    t6 <- paste("The", fn, "function should return a vector of with NaN values when input has inconsistent parameters and not raise an error")
    if (inherits(res3, "try-error")) {
        res$txt <- t6
        return(res)
    }
    x <- 0:1
    names(start.arg) <- paste0(names(start.arg), "_")
    res4 <- try(do.call(fn, c(list(x), start.arg, fix.arg)),      silent = TRUE)
    t8 <- paste("The", fn, "function should raise an error when names are incorrectly named")
    if (!inherits(res4, "try-error")) {
        res$txt <- t8
        return(res)
    }
    return(data.frame(ok = TRUE, txt = ""))
}
