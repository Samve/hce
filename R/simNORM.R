
#' Simulate multivariate normal distribution using compound symmetry covariance structure
#'
#' @param n number of observations.
#' @param mean vector of means.
#' @param sd vector of length 1 of standard deviation.
#' @param rho vector of length 1 of between subject correlation.
#'
#' @return a matrix with `n` rows and `length(mean)` columns.
#' @export
#'
#' @examples
#' MN <- simNORM(n = 1000, mean = c(1, 2), sd = 2, rho = -0.5)
#' colMeans(MN)
#' var(MN)
simNORM <- function(n, mean, sd = 1, rho = 0)
{
  m <- length(mean)
  sd <- sd[1]
  rho <- rho[1]
  
  stopifnot("Correlation rho should be <= 1" = abs(rho) <= 1)
  
  sigma <- matrix(data = sd^2*rep(rho, times = m^2), nrow = m)
  diag(sigma) <- sd^2
  
  Q <- eigen(sigma, symmetric = TRUE)
  
  stopifnot("The covariance matrix is not positive definite" = any(Q$vectors < 0))
  
  sqrtQ <-  Q$vectors %*% diag(x = sqrt(Q$values), nrow = length(Q$values)) %*% t(Q$vectors) 
  X <- matrix(stats::rnorm(n * ncol(sigma)), nrow = n, byrow = TRUE) %*%  sqrtQ
  X <- sweep(x = X, MARGIN = 2, STATS = mean, FUN = "+")   
  colnames(X) <- names(mean)
  X
}
