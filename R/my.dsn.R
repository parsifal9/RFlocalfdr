#' density of skew-normal
#'
#' This function allows you to express your love of cats.
#' @param x param
#' @param xi param
#' @param omega param
#' @param lambda param
#' @keywords cats
#' @export
#' @examples
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
