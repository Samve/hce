#' Simulate random numbers from a Weibull distribution with gamma frailty
#'
#' @param n number of observations.
#' @param shape shape parameter of the Weibull distribution.
#' @param rate rate parameter of the Weibull distribution.
#' @param theta Variance parameter of the gamma frailty distribution, with a default value of 1. A value of 0 corresponds to no frailty, that is, a standard Weibull distribution. 
#' @details
#' The survival function is given by \deqn{S(t|\lambda, \alpha, \gamma) = e^{-\gamma\cdot\lambda\cdot t^\alpha},\ \ \alpha>0,\,\lambda>0,} where \eqn{\alpha} is the shape parameter, \eqn{\lambda} is the rate parameter, and \eqn{\gamma} is the frailty generated from a gamma distribution with \eqn{shape = \frac{1}{\theta}} 
#' and \eqn{scale=\theta} which results in a frailty variance of \eqn{\theta}. When \eqn{\theta = 0}, the distribution reduces to a standard Weibull distribution without frailty.
#' 
#' The marginal (unconditional) survival function of the Weibull distribution with gamma frailty can be expressed as \deqn{S(t|\lambda, \alpha, \theta) = \left(1 + \theta\cdot\lambda\cdot t^\alpha\right)^{-\frac{1}{\theta}},\ \ \alpha>0,\,\lambda>0,\,\theta\geq 0.}
#' @return a vector of random numbers of length `n`.
#' @export
#' @md
#' @references Wienke A. ``Frailty Models in Survival Analysis." Chapman and Hall/CRC (2010).
#' @seealso [hce::simTTE()]  for simulation of a two-event time-to-event endpoint from this distribution under the illness-death model.
#' @examples
#' # Simulate 10 random numbers
#' set.seed(123)
#' rweibullGF(n = 10, shape = 2, rate = 0.5, theta = 0.1)
rweibullGF <- function(n, shape, rate, theta = 1) {
  n <- n[1]
  shape <- shape[1]
  rate <- rate[1]
  theta <- theta[1]
  stopifnot("`n` must be a positive integer" = n > 0,
            "`shape` must be a positive number" = shape > 0,
            "`rate` must be a positive number" = rate > 0,
            "`theta` must be a non-negative number" = theta >= 0)
  n <- round(n)
  # When theta = 0, use a regular Weibull distribution without frailty. When theta is very small, set it to a small positive value to avoid numerical issues.
  ## We use PH (proportional hazards parametization) with rate and shape. The scale parameter in the stats::rweibull function is defined as scale = rate^(-1/shape).
  if(theta == 0)
    return(stats::rweibull(n, shape = shape, scale = rate^(- 1/shape)))
  U <- stats::runif(n)                                      
  G <- stats::rgamma(n, shape = 1/theta, scale = theta)   
  ## Uses PH parametrization for the Weibull survival function
  ## S(t) = exp(-gamma*scale*t^shape) where gamma is the frailty term
  (- log(U) / (G * rate))^(1/shape)                 
}