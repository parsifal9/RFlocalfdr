#' propTrueNullByLocalFDR
#'
#' propTrueNullByLocalFDR
#' @param p  probabilities
#' @keywords cats
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
