#' Simulate random numbers from a Lomax distribution with exponential truncation
#'
#' @param n number of observations.
#' @param rate0 initial hazard rate at time 0. Must be a positive vector of length 1 or `n`.
#' @param rate1 asymptotic hazard rate used as the rate of exponential distribution. Must be a positive vector of length 1 or `n`.
#' @param scale scale parameter of the Lomax distribution (GLL with `shape = 1`). Must be a positive vector of length 1 or `n`.
#' @details
#' The generated random variable \eqn{Z} has hazard function
#' \deqn{
#' h(t) = \lambda_1 +
#' \frac{\lambda_0 - \lambda_1}{1 + t / \tau}, \qquad t \geq 0.
#' }{
#' h(t) = lambda1 + (lambda0 - lambda1) / (1 + t / tau), t >= 0.
#' }
#' where \eqn{\lambda_0 = \texttt{rate0}}, \eqn{\lambda_1 =
#' \texttt{rate1}}, and \eqn{\tau = \texttt{scale}}. The parameters
#' satisfy \eqn{\lambda_0 \geq \lambda_1 > 0} and \eqn{\tau > 0}.
#'
#' At time zero, the hazard is
#' \deqn{
#' h(0) = \lambda_0.
#' }{
#' h(0) = lambda0.
#' }
#' As time tends to infinity, the hazard approaches the asymptotic rate
#' \deqn{
#' \lim_{t \to \infty} h(t) = \lambda_1.
#' }{
#' limit h(t) = lambda1 as t -> infinity.
#' }
#'
#' The cumulative hazard is
#' \deqn{
#' H(t) = \lambda_1 t +
#' (\lambda_0 - \lambda_1)\tau
#' \log\left(1 + \frac{t}{\tau}\right).
#' }{
#' H(t) = lambda1 * t + (lambda0 - lambda1) * tau *
#' log(1 + t / tau).
#' }
#' and the survival function is
#' \deqn{
#' S(t) = \Pr(Z > t) = \exp\{-H(t)\}.
#' }{
#' S(t) = P(Z > t) = exp(-H(t)).
#' }
#'
#' The algorithm generates
#' \deqn{
#' Z = \min(X, Y).
#' }{
#' Z = min(X, Y).
#' }
#' where \eqn{Y \sim \operatorname{Exponential}(\lambda_1)} and \eqn{X}
#' follows a generalized log-logistic distribution with shape parameter
#' 1, rate parameter \eqn{\lambda_0 - \lambda_1}, and
#' \deqn{
#' \theta = \frac{1}{(\lambda_0 - \lambda_1)\tau}.
#' }{
#' theta = 1 / ((lambda0 - lambda1) * tau).
#' }
#' Its survival function is
#' \deqn{
#' S_X(t) = \Pr(X > t)
#' = \left(
#' 1 + \frac{(\lambda_0 - \lambda_1)t}
#' {(\lambda_0 - \lambda_1)\tau}
#' \right)^{-(\lambda_0 - \lambda_1)\tau}
#' = \left(1 + \frac{t}{\tau}\right)^{-(\lambda_0 - \lambda_1)\tau}.
#' }{
#' S_X(t) = P(X > t) =
#' (1 + ((lambda0 - lambda1) * t) /
#' ((lambda0 - lambda1) * tau))^(-(lambda0 - lambda1) * tau) =
#' (1 + t / tau)^(-(lambda0 - lambda1) * tau).
#' }
#' When \eqn{\lambda_0 = \lambda_1}, the Lomax component vanishes and
#' \eqn{Z} is exponential with rate \eqn{\lambda_1}.
#'
#' For independent random variables \eqn{X} and \eqn{Y}, the survival
#' function of their minimum is
#' \deqn{
#' \Pr(Z > t) = \Pr(X > t, Y > t) = S_X(t)S_Y(t).
#' }{
#' P(Z > t) = P(X > t, Y > t) = S_X(t) * S_Y(t).
#' }
#' Therefore, the cumulative hazard of \eqn{Z} is additive:
#' \deqn{
#' H_Z(t) = -\log S_Z(t)
#' = -\log S_X(t) - \log S_Y(t)
#' = H_X(t) + H_Y(t).
#' }{
#' H_Z(t) = -log(S_Z(t)) = H_X(t) + H_Y(t).
#' }
#' Since \eqn{Y} is exponential, \eqn{H_Y(t) = \lambda_1 t}. For the
#' generalized log-logistic component with shape 1,
#' \deqn{
#' H_X(t) =
#' (\lambda_0 - \lambda_1)\tau
#' \log\left(1 + \frac{t}{\tau}\right).
#' }{
#' H_X(t) = (lambda0 - lambda1) * tau *
#' log(1 + t / tau).
#' }
#' Adding these cumulative hazards gives \eqn{H(t)}. Differentiating
#' \eqn{H(t)} gives the stated hazard function. 
#' 
#' @return a vector of random numbers of length `n`. When vectors of length `n` are provided, each element of the output is generated using the corresponding elements of `rate0`, `rate1`, and `scale`. 
#' @export
#' @md
#' @references Lomax KS. "Business failures: Another example of the analysis
#' of failure data." Journal of the American Statistical Association 49.268
#' (1954): 847-852. <doi:10.1080/01621459.1954.10501239>.
#' @seealso [hce::GLL()] for the Generalized Log-Logistic, including Lomax distribution.
#' @examples
#' rlomaxET(10, rate0 = 0.12, rate1 = 0.1, scale = 1)
rlomaxET <- function(n, rate0, rate1, scale = 1){
  # Lomax random variable with exponential truncation
  # Random variable with hazard:
  # h(t) = rate1 + (rate0 - rate1) / (1 + t / scale)
  n <- as.integer(n)[1]
  stopifnot("The sample size `n` must be a postive integer." = length(n) == 1L && !is.na(n) && n > 0L,
            "The initial event rate `rate0` must be a positive vector." = is.numeric(rate0) && all(rate0 > 0),
            "The final event rate `rate1` must be a positive vector." =  is.numeric(rate1) && all(rate1 > 0),
            "The scale parameter `scale` must be a positive vector." =  is.numeric(scale) && all(scale > 0),
            "The length of the vector `rate0` must be equal to 1 or `n`" = length(rate0) %in% c(1, n),
            "The length of the vector `rate1` must be equal to 1 or `n`" = length(rate1) %in% c(1, n),
            "The length of the vector `scale` must be equal to 1 or `n`" = length(scale) %in% c(1, n))
  
  rate0 <- rep_len(rate0, n)
  rate1 <- rep_len(rate1, n)
  scale <- rep_len(scale, n)
  
  stopifnot(
    "The initial rate `rate0` must be greater than or equal to the asymptotic rate `rate1`." = 
      all(rate0 >= rate1)
  )
  
  rate <- rate0 - rate1
  Z <- numeric(n)
  # Exponential component
  Y <- stats::rexp(n, rate = rate1)
  
  # Cases where rate0 == rate1: only the exponential component remains
  equal_rate <- rate == 0
  Z[equal_rate] <- Y[equal_rate]
  
  # Cases where rate0 > rate1: combine the two hazards
  unequal_rate <- rate > 0
  
  if(any(unequal_rate)){
    theta <- 1/(rate[unequal_rate] * scale[unequal_rate])
    Z[unequal_rate] <- pmin(rGLL(sum(unequal_rate), rate = rate[unequal_rate], shape = 1, theta = theta ),
                            Y[unequal_rate])
  }
  Z
}
