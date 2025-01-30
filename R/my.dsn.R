#' my.dsn
#'
#' density of skew-normal using the appromimation of Ashour, Samir K. and Abdel-hameed, Mahmood A.
#'
#' See \url{https://en.wikipedia.org/wiki/Skew_normal_distribution} for discussion of the skew-normal.
#' Using the appromimation of Ashour, Samir K. and Abdel-hameed, Mahmood A.
#' "Approximate Skew Normal Distribution", Journal of Advanced Research, 2010, 1:4.
#' It accepts the parameters xi, omega, lambda (Ashour et. al. 2010). Other formulations may use
#' different parameterizations. The sn (skew-normal) package incluse the extended skew-normal (ESN) distribution. For the SN the
#' tau parameter is 0.
#'
#' RFlocalfdr also uses wrappers around the functions dsn, qsn and psn from the package "sn"
#' https://cran.r-project.org/web/packages/sn/.
#' This is due to the fact that 
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
#' \dontrun{
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


my.dsn <- function(x, xi=0, omega=1, lambda=1) {
    if(abs(lambda) > 5) {
        warning("The approximation may be unreliable for |lambda| > 5")
    }
    
    # For negative lambda, reflect x around xi and use |lambda|
    if(lambda < 0) {
        return(my.dsn(2*xi - x, xi=xi, omega=omega, lambda=abs(lambda)))
    }
    
    z <- (x - xi)/omega
    aa <- (1/sqrt(2*pi))*exp(-z^2/2)
    bb <- 1/3*lambda^3*z^3
    cc <- 3*lambda^2*z^2
    
    # Initialize output vector
    hx <- numeric(length(x))
    
    # Fill in values for each region
    idx_1 <- z < -3/lambda
    idx_2 <- z >= -3/lambda & z < -1/lambda
    idx_3 <- z >= -1/lambda & z < 1/lambda
    idx_4 <- z >= 1/lambda & z < 3/lambda
    idx_5 <- z >= 3/lambda
    
    hx[idx_1] <- 0  # Changed from .Machine$double.xmin to 0
    hx[idx_2] <- (1/8)*aa[idx_2]*(9*lambda*z[idx_2] + cc[idx_2] + bb[idx_2] + 9)
    hx[idx_3] <- (1/4)*aa[idx_3]*(3*lambda*z[idx_3] + bb[idx_3] + 4)
    hx[idx_4] <- (1/8)*aa[idx_4]*(9*lambda*z[idx_4] - cc[idx_4] + bb[idx_4] + 7)
    hx[idx_5] <- sqrt(2/pi)*exp(-z[idx_5]^2/2)  # Fixed this line
    
    (1/omega)*hx
}

## my.dsn<-function(x,xi=0, omega = 1,lambda=1){
##     z <- (x - xi)/omega
##     aa<- (1/sqrt(2*pi))*exp(-z^2/2)
##     bb<- 1/3*lambda^3*z^3
##     cc <- 3*lambda^2*z^2
##     hx<-
##         ifelse( (z < -3/lambda),0,
##         ifelse((z < -1/lambda),(1/8)*aa*(9*lambda*z + cc+ bb+9),
##         ifelse( (z  < 1/lambda), (1/4)*aa*(3*lambda*z + bb + 4),
##         ifelse( (z < 3/lambda), (1/8)*aa*(9*lambda*z -cc +bb+7),sqrt(2/pi)*exp(-z^2/2)
##         ))))
##     (1/omega)*hx
## }

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
     aa<- (
        (names(match.call())[3]=="xi") &
        (names(match.call())[4]=="omega") &
        (names(match.call())[5]=="alpha") &
        (names(match.call())[6]=="tau")
    )
    if (!(aa)){
        stop("error! names are incorrect")
    }

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

    if (any(is.na(p))){
        return(rep(NaN,length(p)))
        }

    if (any(is.nan(p))){
        return(rep(NaN,length(p)))
        }
    
    if (any(p <0) | any(p>1)){
        return(rep(NaN,length(p)))
    }

    
    
    sn::qsn(x=p, xi = xi, omega =omega, alpha=alpha ,tau = tau,dp=NULL)

}
