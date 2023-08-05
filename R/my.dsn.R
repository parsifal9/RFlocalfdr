#' my.dsn
#'
#' density of skew-normal using the appromimation of Ashour, Samir K. and Abdel-hameed, Mahmood A.
#'
#' See \url{https://en.wikipedia.org/wiki/Skew_normal_distribution} for discussion of the skew-normal.
#' Using the appromimation of Ashour, Samir K. and Abdel-hameed, Mahmood A.
#' "Approximate Skew Normal Distribution", Journal of Advanced Research, 2010, 1:4.
#' It accepts the parameters xi, omega, lambda (Ashour et. al. 2010). Other foumulations may use
#' different parameterizations. The sn (skew-normal) package incluse the extended skew-normal (ESN) distribution. For the SN the
#' tau parameter is 0.
#'
#' fitdistrplus::fitdist(imp, "sn", start = list(xi = mean(imp)...
#' returns warnings such as
#' The dsn function should return a zero-length vector when input has length zero and not raise an error
#' The psn function should have its first argument named: q as in base R
#' These wrappers ensure conformity with the expectations of fitdistrplus::fitdist
#' 
#' @param x vector of quantiles. Missing values (‘NA’'s) and ‘Inf’'s are   allowed.
#' @param p vector of probabilities. Missing values (‘NA’'s) are allowed
#' @param q vector of quantiles
#' @param xi vector of location parameters.
#' @param alpha vector of slant parameter(s)
#' @param omega vector of scale parameters; must be positive
#' @param tau = 0
#' @param lambda param
#' @param ... arguments passed to sn
#' @keywords skew normal
#' @export my.dsn
#' @export dsn
#' @export psn
#' @export qsn
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

my.dsn<-function(x,xi=0, omega = 1,lambda=1){
    z <- (x - xi)/omega
    aa<- (1/sqrt(2*pi))*exp(-z^2/2)
    bb<- 1/3*lambda^3*z^3
    cc <- 3*lambda^2*z^2
    hx<-
        ifelse( (z < -3/lambda),0,
        ifelse((z < -1/lambda),(1/8)*aa*(9*lambda*z + cc+ bb+9),
        ifelse( (z  < 1/lambda), (1/4)*aa*(3*lambda*z + bb + 4),
        ifelse( (z < 3/lambda), (1/8)*aa*(9*lambda*z -cc +bb+7),sqrt(2/pi)*exp(-z^2/2)
        ))))
    (1/omega)*hx
}

#' @rdname my.dsn 
dsn<-function(x, xi = 0, omega = 1, alpha = 0, tau = 0){
    if (length(x)==0){
        return( x)
    }
    aa<- (
        (names(match.call())[3]=="xi") &
        (names(match.call())[4]=="omega") &
        (names(match.call())[5]=="alpha") &
        (names(match.call())[6]=="tau")
    )
    if (!(aa)){
        stop("error! names are incorrect")
    }
    res <- sn::dsn(x=x, xi = xi, omega =omega, alpha=alpha ,tau = tau,dp=NULL,log=FALSE)
    res
}

#' @rdname my.dsn
psn <- function(q, xi = -Inf, omega =1, alpha=0 ,tau = 0,...) {
    x <- q

    aa<- (
        (names(match.call())[3]=="xi") &
        (names(match.call())[4]=="omega") &
        (names(match.call())[5]=="alpha") &
        (names(match.call())[6]=="tau")
    )
    
    if (!(aa)){
        stop("error! names are incorrect")
    }
    
    sn::psn(x, xi = xi, omega =omega, alpha=alpha ,tau = tau,dp=NULL)
}

#' @rdname my.dsn
qsn <- function(p, xi = Inf, omega =1, alpha=0 ,tau = 0,...) {
    inconsistent <- FALSE
    if (xi >0){
        inconsistent <- TRUE
    }
    if(omega < 0){
        inconsistent <- TRUE
    }

    if (inconsistent){
        return(rep(NaN,length(p)))
    }
    
    if (length(p)==0){
        return( p)
    }

    if (any(p <0) |  any(is.nan(p))| any(p>1)){
        return(rep(NaN,length(p)))
    }
    
    sn::qsn(x=p, xi = xi, omega =omega, alpha=alpha ,tau = tau,dp=NULL)

}


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
