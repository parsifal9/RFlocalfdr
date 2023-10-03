#' propTrueNullByLocalFDR
#'
#' Estimate proportion of NULL p-values. Based on  .propTrueNullByLocalFDR in  limma: Linear Models for Microarray Data
#' written by Belinda Phipson and Gordon Smyth
#' @param p  probabilities
#' @keywords cats
#' @return An estimate of the proportion of null p-values by the local fdr
#' @export
propTrueNullByLocalFDR <- function(p) 
{
    n <- length(p)
    i <- n:1L
    p <- sort(p, decreasing = TRUE)
    q <- pmin(n/i * p, 1)
    n1 <- n + 1L
    sum(i * q)/n/n1 * 2
 }
