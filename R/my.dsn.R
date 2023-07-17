#' density of skew-normal using the appromimation of Ashour, Samir K. and Abdel-hameed, Mahmood A.
#'
#' See \url{https://en.wikipedia.org/wiki/Skew_normal_distribution} for discussion of the skew-normal.
#' Using the appromimation of Ashour, Samir K. and Abdel-hameed, Mahmood A.
#' "Approximate Skew Normal Distribution", Journal of Advanced Research, 2010, 1:4.
#' It accepts the parameters xi, omega, lambda (Ashour et. al. 2010). Other foumulations may use
#' different parameterizations. The sn (skew-normal) package incluse the extended skew-normal (ESN) distribution. For the SN the
#' tau parameter is 0.
#' 
#' @param x data
#' @param xi param
#' @param omega param
#' @param lambda param
#' @keywords skew normal
#' @export
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



