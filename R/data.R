#' 20000 importance values
#'
#' A dataset containing  20000 importance values
#'
#' @format A vector varaible importances with 20000 values
#' \describe{
#'   \item{imp1}{importances }
#' }
#' @examples
#' \donttest{
#' require(ranger)
#' inv.logit <-function (x) {
#'     plogis(x)}
#' 
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
#' set.seed(19)
#' X<- make_data(nVars,nSamples)
#'
#' w <- rep(0, times = nVars)
#' w[101] <- 1
#' w[102] <- 1/sqrt(2)
#' w[103] <- 1/sqrt(4)
#' w[104] <- 1/sqrt(8)
#' w[105] <- 1/sqrt(16)
#' y <- make_response(X, w)
#'
#'
#' colnames(X) <- c(make.names(1:20000))
#' set.seed(19)
#' rf1<-ranger::ranger(y=y,x=X, num.trees = 2000,importance="impurity")
#' table(y,predict(rf1,data=X)$predictions)
#' #OOB prediction error:             41.30 % 
#' table(y,predict(rf1,data=X)$predictions)
#'
#' t2 <-count_variables(rf1)
#' head(t2)
#' dim(t2)
#'
#' imp<-rf1$variable.importance
#' imp<-log(imp)
#' plot(density((imp)))
#' hist(imp,col=6,lwd=2,breaks=100,main="histogram of importances")
#'
#' res.temp <- determine_cutoff(imp, t2 ,cutoff=c(0,1,2,3),plot=c(0,1,2,3),Q=0.75,try.counter=1)
#' plot(c(0,1,2,3),res.temp[,3])
#' imp<-imp[t2 > 1]
#' qq <- plotQ(imp,debug.flag = 0)
#' ppp<-run.it.importances(qq,imp,debug=0)
#' aa<-significant.genes(ppp,imp,cutoff=0.2,debug.flag=0,do.plot=2, use_95_q=TRUE)
#' length(aa$probabilities)
#' names(aa$probabilities)
#' #' #[1] "X101"   "X102"   "X103"   "X104"   "X105"   "X2994"  "X9365"  "X10718"
#' # [9] "X13371" "X15517" "X16460"
#'
#' counts<-t2
#' imp20000 <- list(imp,counts)
#' names(imp20000) <-c("importances","counts")
#' }
#' 
"imp20000"
