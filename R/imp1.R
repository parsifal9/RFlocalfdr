#' 2000 importance values
#'
#' A dataset containing  2000 importance values
#'
#' @format A vector with 2000 values
#' \describe{
#'   \item{imp1}{importaces }
#' }
#' @source \url{http://www.diamondse.info/}
#' @examples
#' library(ranger)
#' inv.logit <-function (x) {
#'     plogis(x)}
#' make_data <- function(nVars, nSamples) {
#'     as.matrix(sapply(1:nVars, function(t){sample(0:2, nSamples, replace=TRUE)}))
#' }
#' 
#' make_cont_response <- function(X, w) {
#'     (X-1) %*% w 
#' }
#' 
#' make_response <- function(X, w) {
#'    as.factor(inv.logit((X-1) %*% w * 2 ) > runif(nrow(X)))
#' } 
#' 
#' 
#' nVars <-   20000
#' nSamples <- 1000
#' 
#'set.seed(19)
#'X<- make_data(nVars,nSamples)
#'
#'w <- rep(0, times = nVars)
#'w[101] <- 1
#'w[102] <- 1/sqrt(2)
#'w[103] <- 1/sqrt(4)
#'w[104] <- 1/sqrt(8)
#'w[105] <- 1/sqrt(16)
#'X.temp<- make_data(nVars,nSamples)
#'y <- make_response(X.temp, w)
#'
#'
#'colnames(X) <- c(make.names(1:20000))
#'library(ranger)
#'rf1<-ranger(y=y,x=X, num.trees = 20000,importance="impurity")
#'table(y,predict(rf1,data=X)$predictions)
#'#OOB prediction error:             51.70 % 
#'table(y,predict(rf1,data=X)$predictions)
#'
#'system.time(t2<- table(factor(aa1[aa2 != 0],levels=1:10000))) # 6.620   
#'head(t2)
#'dim(t2)
#'
#'imp<-rf1$variable.importance
#'
#'imp<-log(imp)
#'plot(density((imp)))
#'hist(imp,col=6,lwd=2,breaks=100,main="histogram of importances")
