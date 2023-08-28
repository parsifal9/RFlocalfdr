#' 2000 importance values
#'
#' A dataset containing  2000 importance values
#'
#' @format A vector with 2000 values
#' \describe{
#'   \item{imp1}{importaces }
#' }
#' @examples
#' \dontrun{
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
#' set.seed(19)
#' X<- make_data(nVars,nSamples)
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
#'rf1<-ranger::ranger(y=y,x=X, num.trees = 20000,importance="impurity")
#'table(y,predict(rf1,data=X)$predictions)
#'#OOB prediction error:             51.70 % 
#'table(y,predict(rf1,data=X)$predictions)
#'
#'t2 <-count_variables(rf1)
#'head(t2)
#'dim(t2)
#'
#'imp<-rf1$variable.importance
#'imp<-log(imp)
#'plot(density((imp)))
#'hist(imp,col=6,lwd=2,breaks=100,main="histogram of importances")
#' }
"imp1"


## #' ch22 importance values
## #'
## #' A dataset containing  1103547 importance values, and a table of variables used in splits
## #'
## #' @format A list
## #' \describe{
## #'   \item{imp}{importaces }
## #'   \item{C}{table of counts}
## #' }
## #' @source "A Global Reference for Human Genetic Variation", Auton et al., Nature, 2015, 526:7571 pp 68--74
## #' @examples
## #' \dontrun{
## #' library(ranger)
## #' system.time(fit.ranger.7 <- ranger(dependent.variable.name= "V1", data = aa2,
## #'                                 importance = "impurity",
## #'                                  num.threads=20,num.trees = 100000,
## #'                                  seed=123))
## #' }
## #' #Ranger result
## #' #Call:
## #' #ranger(dependent.variable.name = "V1", data = aa2, importance = "impurity", 
## #' #                              num.threads = 20, num.trees = 1e+05, seed = 123) 
## #' #Type:                             Classification 
## #' #Number of trees:                  1e+05 
## #' #Sample size:                      2504 
## #' #Number of independent variables:  1103547 
## #' #Mtry:                             1050 
## #' #Target node size:                 1 
## #' #Variable importance mode:         impurity 
## #' #Splitrule:                        gini 
## #' #OOB prediction error:             4.27 %
## "ch22"

## #' Effects of cigarette smoke on the human airway epithelial cell transcriptome
## #' 
## #' A dataset containing normalized transcript measurements for 51 subjects and 22283 transcripts.
## #' See Spira et al (2004). "Gene Expression Profiling of Human Lung Tissue from Smokers with Severe Emphysema",
## #' Am J Respir Cell Mol Biol.
## #'
## #' @format A list with rma (the transcript data) and y (the class labels):
## #' \describe{
## #'   \item{rma}{ 51 by  22283, log2 real values }
## #'   \item{y}{a character vector, "smoking" and "never-smoked"}
## #'   ...
## #' }
## #' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE994}
## "smoking"
