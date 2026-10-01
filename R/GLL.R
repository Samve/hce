#' The Generalized Log-Logistic Distribution
#'
#' Density, distribution function, quantile function, hazard, cumulative hazard, survival odds,
#' and random generation for the generalized log-logistic distribution.
#'
#' @param x,q non-negative values at which to evaluate the density,
#'   distribution, hazard, cumulative hazard, or survival odds.
#' @param p probabilities in \eqn{[0, 1]}. When `log.p = TRUE`, `p` contains
#'   log-probabilities in \eqn{[-\infty, 0]}.
#' @param n number of observations. Must be a single non-negative integer.
#' @param rate positive rate parameter of the distribution.
#' @param shape positive shape parameter of the distribution.
#' @param theta non-negative additional parameter of the distribution, with
#'   default value `0`.
#' @param log logical; if `TRUE`, the density is returned on the log scale.
#' @param log.p logical; if `TRUE`, probabilities are supplied or returned
#'   on the log scale.
#' @param lower.tail logical; if `TRUE` (default), probabilities are
#'   \eqn{P[X \leq x]}; otherwise, they are \eqn{P[X > x]}.
#' @name GLL
#'
#' @details
#' The generalized log-logistic distribution with rate parameter \eqn{\lambda > 0},
#' shape parameter \eqn{\alpha > 0}, and parameter \eqn{\theta \ge 0} has survival
#' function
#' \deqn{
#' S(t \mid \lambda, \alpha, \theta) =
#' \left(1 + \theta \lambda t^\alpha\right)^{-1 / \theta},
#' \qquad \theta > 0.
#' }
#'
#' The corresponding cumulative distribution function is
#' \deqn{
#' F(t \mid \lambda, \alpha, \theta) =
#' 1 - \left(1 + \theta \lambda t^\alpha\right)^{-1 / \theta}.
#' }
#'
#' For \eqn{\theta > 0}, the quantile function for a lower-tail
#' probability \eqn{p = F(t)} is
#' \deqn{
#' Q(p \mid \lambda, \alpha, \theta) =
#' \left(
#' \frac{(1 - p)^{-\theta} - 1}{\theta \lambda}
#' \right)^{1 / \alpha},
#' \qquad 0 \leq p \leq 1.
#' }
#'
#' When `lower.tail = FALSE`, `p` is an upper-tail probability
#' \eqn{p = S(t)}, and the quantile is
#' \deqn{
#' Q(p \mid \lambda, \alpha, \theta) =
#' \left(
#' \frac{p^{-\theta} - 1}{\theta \lambda}
#' \right)^{1 / \alpha},
#' \qquad 0 \leq p \leq 1.
#' }
#'
#' When `log.p = TRUE`, `p` is the logarithm of the selected
#' lower- or upper-tail probability.
#'
#' The cumulative hazard function is
#' \deqn{
#' H(t \mid \lambda, \alpha, \theta) =
#' \frac{\log\left(1 + \theta \lambda t^\alpha\right)}{\theta}.
#' }
#'
#' The density function is
#' \deqn{
#' f(t \mid \lambda, \alpha, \theta) =
#' \lambda \alpha t^{\alpha - 1}
#' \left(1 + \theta \lambda t^\alpha\right)^{-1 / \theta - 1}.
#' }
#'
#' For `dGLL(..., log = TRUE)`, the logarithm of the density is returned. For
#' \eqn{\theta > 0},
#' \deqn{
#' \log f(t \mid \lambda, \alpha, \theta) =
#' \log(\lambda) + \log(\alpha) + (\alpha - 1)\log(t)
#' - \left(\frac{1}{\theta} + 1\right)\log\left(1 + \theta \lambda t^\alpha\right).
#' }
#'
#' The hazard function is
#' \deqn{
#' h(t \mid \lambda, \alpha, \theta) =
#' \frac{\lambda \alpha t^{\alpha - 1}}
#' {1 + \theta \lambda t^\alpha}.
#' }
#'
#' Random generation is performed by inverse transform sampling. If
#' \eqn{U \sim \mathrm{Uniform}(0,1)}, then
#' \deqn{
#' T = Q(U \mid \lambda, \alpha, \theta)
#' }
#' has the generalized log-logistic distribution.
#'
#' As \eqn{\theta \to 0}, this reduces to the Weibull distribution with survival
#' function
#' \deqn{
#' S(t \mid \lambda, \alpha) = \exp\left(-\lambda t^\alpha\right).
#' }
#' The lower-tail quantile function is
#' \deqn{
#' Q(p \mid \lambda, \alpha) =
#' \left(
#' \frac{-\log(1 - p)}{\lambda}
#' \right)^{1 / \alpha},
#' \qquad 0 \leq p \leq 1.
#' }
#'
#' For an upper-tail probability \eqn{p = S(t)}, the corresponding quantile is
#' \deqn{
#' Q(p \mid \lambda, \alpha) =
#' \left(
#' \frac{-\log(p)}{\lambda}
#' \right)^{1 / \alpha}.
#' }
#'
#' When `log.p = TRUE`, `p` is the logarithm of the selected lower- or
#' upper-tail probability.
#'
#' In this case, `dGLL(..., log = TRUE)` returns the corresponding Weibull log-density.
#' 
#' This corresponds to the Weibull distribution in [stats::dweibull()],
#' [stats::pweibull()], [stats::qweibull()], and [stats::rweibull()] with
#' `shape = shape` and `scale = rate^(-1 / shape)`.
#'
#' When \eqn{\theta = 1}, the survival function becomes
#' \deqn{
#' S(t \mid \lambda, \alpha, 1) = \left(1 + \lambda t^\alpha\right)^{-1},
#' }
#' corresponding to the log-logistic distribution.
#' 
#' The survival odds are defined as
#' \deqn{
#' SO(t \mid \lambda, \alpha, \theta)
#' = \frac{S(t \mid \lambda, \alpha, \theta)}
#' {1 - S(t \mid \lambda, \alpha, \theta)}.
#' }
#'
#' For \eqn{\theta > 0},
#' \deqn{
#' SO(t \mid \lambda, \alpha, \theta)
#' =
#' \left[
#' \left(1 + \theta \lambda t^\alpha\right)^{1/\theta}
#' - 1
#' \right]^{-1}.
#' }
#'
#' When \eqn{\theta = 0}, the survival odds are given by the
#' Weibull limit
#' \deqn{
#' SO(t \mid \lambda, \alpha, 0)
#' =
#' \left[\exp\left(\lambda t^\alpha\right)-1\right]^{-1}.
#' }
#'
#' @return
#' For `dGLL()`, `pGLL()`, `qGLL()`, `hGLL()`, `HGLL()`, and `oGLL()`, a numeric
#' vector with the same length as the main input (`x`, `q`, or `p`). For
#' `rGLL()`, a numeric vector of length `n`. Scalar values of `rate`, `shape`,
#' and `theta` are recycled to match the length of the main input. Non-scalar
#' parameter vectors must have the same length as the main input. When
#' vector-valued parameters are supplied, each output element is calculated
#' using the corresponding parameter values.
#' 
#' @md
#'
#' @references
#' Wienke A. *Frailty Models in Survival Analysis*. Chapman and Hall/CRC (2010).
#'
#' @seealso
#' [hce::rweibullGF()] for simulation of a Weibull distribution with gamma frailty.
#'
#' @examples
#' ## Example 1: Compare the CDF of GLL with different values of `theta`
#' shape <- 0.85
#' rate <- 0.1
#' x <- seq(0, 3, length.out = 100)
#' y0 <- pGLL(x, rate = rate, shape = shape, theta = 0)
#' y1 <- pGLL(x, rate = rate, shape = shape, theta = 1)
#' y2 <- pGLL(x, rate = rate, shape = shape, theta = 4)
#' plot(x, y0, type = "l", ylab = "Cumulative distribution function", ylim = c(0, .3), 
#' main = "Generalized Log-Logistic Distribution CDF")
#' lines(x, y1, col = "red", lty = 4, lwd = 2)
#' lines(x, y2, col = "blue", lty = 5, lwd = 2)
#' legend("topleft", legend = c("theta = 0", "theta = 1", "theta = 4"), 
#' col = c("black", "red", "blue"), lty = c(1, 4, 5), lwd = c(1, 2, 2))
NULL

#' @export
#' @rdname GLL
dGLL <- function(x, rate, shape = 1, theta = 0, log = FALSE) {
  n <- length(x)
  log <- as.logical(log[1])
  
  stopifnot(
    "`shape` must be a positive vector of length 1 or length of `x`." =
      is.numeric(shape) && all(shape > 0) && length(shape) %in% c(1, n),
    "`rate` must be a positive vector of length 1 or length of `x`." =
      is.numeric(rate) && all(rate > 0) && length(rate) %in% c(1, n),
    "`theta` must be a non-negative vector of length 1 or length of `x`." =
      is.numeric(theta) && all(theta >= 0) && length(theta) %in% c(1, n),
    "`x` quantiles must be a non-negative vector." =
      is.numeric(x) && all(x >= 0)
  )
  
  rate <- rep_len(rate, length.out = n)
  shape <- rep_len(shape, length.out = n)
  theta <- rep_len(theta, length.out = n)
  
  out <- numeric(n)
  i0 <- theta == 0
  i1 <- theta > 0
  
  if (any(i0)) {
    out[i0] <- stats::dweibull(
      x[i0],
      shape = shape[i0],
      scale = rate[i0]^(-1 / shape[i0]),
      log = log
    )
  }
  
  if (any(i1)) {
    if (log) {
      out[i1] <- log(rate[i1]) + log(shape[i1]) -
        (1 / theta[i1] + 1) * log1p(theta[i1] * rate[i1] * x[i1]^shape[i1])
      
      iz <- i1 & x == 0
      ip <- i1 & x > 0
      
      out[ip] <- out[ip] + (shape[ip] - 1) * log(x[ip])
      
      if (any(iz)) {
        out[iz] <- ifelse(
          shape[iz] < 1,
          Inf,
          ifelse(shape[iz] == 1, log(rate[iz]) + log(shape[iz]), -Inf)
        )
      }
    } else {
      out[i1] <- rate[i1] * shape[i1] * x[i1]^(shape[i1] - 1) *
        (1 + theta[i1] * rate[i1] * x[i1]^shape[i1])^(-1 / theta[i1] - 1)
    }
  }
  
  out
}

#' @export
#' @rdname GLL
pGLL <- function(q, rate, shape = 1, theta = 0, lower.tail = TRUE, log.p = FALSE) {
  n <- length(q)
  lower_tail <- as.logical(lower.tail[1])
  log_p <- as.logical(log.p[1])
  
  stopifnot(
    "`shape` must be a positive vector of length 1 or length of `q`." =
      is.numeric(shape) && all(shape > 0) && length(shape) %in% c(1, n),
    "`rate` must be a positive vector of length 1 or length of `q`." =
      is.numeric(rate) && all(rate > 0) && length(rate) %in% c(1, n),
    "`theta` must be a non-negative vector of length 1 or length of `q`." =
      is.numeric(theta) && all(theta >= 0) && length(theta) %in% c(1, n),
    "`q` quantiles must be a non-negative vector." =
      is.numeric(q) && all(q >= 0)
  )
  
  rate <- rep_len(rate, length.out = n)
  shape <- rep_len(shape, length.out = n)
  theta <- rep_len(theta, length.out = n)
  
  out <- numeric(n)
  i0 <- theta == 0
  i1 <- theta > 0
  
  if (any(i0)) {
    out[i0] <- stats::pweibull(
      q[i0],
      shape = shape[i0],
      scale = rate[i0]^(-1 / shape[i0]),
      lower.tail = lower_tail,
      log.p = log_p
    )
  }
  
  if (any(i1)) {
    # Cumulative Hazard
    H <- log1p(
      theta[i1] * rate[i1] * q[i1]^shape[i1]
    ) / theta[i1]
    
    if (log_p) {
      if (lower_tail) {
        # log(1 - exp(-H)), evaluated stably
        out[i1] <- log(-expm1(-H))
      } else {
        # log(S) = -H
        out[i1] <- -H
      }
    } else {
      if (lower_tail) {
        # 1 - exp(-H), evaluated stably
        out[i1] <- -expm1(-H)
      } else {
        out[i1] <- exp(-H)
      }
    }
  }
  
  out
}

#' @export
#' @rdname GLL
qGLL <- function(p, rate, shape = 1, theta = 0, lower.tail = TRUE, log.p = FALSE) {
  n <- length(p)
  lower_tail <- as.logical(lower.tail[1])
  log_p <- as.logical(log.p[1])
  
  stopifnot(
    "`shape` must be a positive vector of length 1 or length of `p`." =
      is.numeric(shape) &&
      all(shape > 0) &&
      length(shape) %in% c(1, n),
    
    "`rate` must be a positive vector of length 1 or length of `p`." =
      is.numeric(rate) &&
      all(rate > 0) &&
      length(rate) %in% c(1, n),
    
    "`theta` must be a non-negative vector of length 1 or length of `p`." =
      is.numeric(theta) &&
      all(theta >= 0) &&
      length(theta) %in% c(1, n),
    
    "`p` must be numeric; values must be in [0, 1] when `log.p = FALSE` and <= 0 when `log.p = TRUE`." =
      if (log_p) {
        is.numeric(p) && all(p <= 0)
      } else {
        is.numeric(p) && all(p >= 0) && all(p <= 1)
      }
  )
  
  rate <- rep_len(rate, length.out = n)
  shape <- rep_len(shape, length.out = n)
  theta <- rep_len(theta, length.out = n)
  
  out <- numeric(n)
  i0 <- theta == 0
  i1 <- theta > 0
  
  if (any(i0)) {
    out[i0] <- stats::qweibull(
      p[i0],
      shape = shape[i0],
      scale = rate[i0]^(-1 / shape[i0]),
      lower.tail = lower_tail,
      log.p = log_p
    )
  }
  
  if (any(i1)) {
    # Calculate the cumulative hazard H = -log(S).
    #
    # For lower-tail input:
    #   p = F(q)
    #
    # For upper-tail input:
    #   p = S(q)
    H <- if (lower_tail) {
      if (log_p) {
        # p is log(F); log(S) = log(1 - F)
        -log(-expm1(p[i1]))
      } else {
        # p is F; log1p(-p) is stable when p is close to zero
        -log1p(-p[i1])
      }
    } else {
      if (log_p) {
        # p is log(S)
        -p[i1]
      } else {
        # p is S
        -log(p[i1])
      }
    }
    
    # q = [(expm1(theta * H)) / (theta * rate)]^(1 / shape)
    z <- expm1(theta[i1] * H)
    
    out[i1] <- (z / (theta[i1] * rate[i1]))^(1 / shape[i1])
  }
  
  out
}

#' @export
#' @rdname GLL
rGLL <- function(n, rate, shape = 1, theta = 0) {
  stopifnot(
    "`n` must be a single non-negative integer." =
      is.numeric(n) && length(n) == 1 && n >= 0 && n == as.integer(n),
    "`shape` must be a positive vector of length 1 or length `n`." =
      is.numeric(shape) && all(shape > 0) && length(shape) %in% c(1, n),
    "`rate` must be a positive vector of length 1 or length `n`." =
      is.numeric(rate) && all(rate > 0) && length(rate) %in% c(1, n),
    "`theta` must be a non-negative vector of length 1 or length `n`." =
      is.numeric(theta) && all(theta >= 0) && length(theta) %in% c(1, n)
  )
  
  qGLL(
    p = stats::runif(n),
    rate = rate,
    shape = shape,
    theta = theta
  )
}

#' @export
#' @rdname GLL
hGLL <- function(x, rate, shape = 1, theta = 0) {
  n <- length(x)
  stopifnot(
    "`shape` must be a positive vector of length 1 or length of `x`." =
      is.numeric(shape) && all(shape > 0) && length(shape) %in% c(1, n),
    "`rate` must be a positive vector of length 1 or length of `x`." =
      is.numeric(rate) && all(rate > 0) && length(rate) %in% c(1, n),
    "`theta` must be a non-negative vector of length 1 or length of `x`." =
      is.numeric(theta) && all(theta >= 0) && length(theta) %in% c(1, n),
    "`x` quantiles must be a non-negative vector." =
      is.numeric(x) && all(x >= 0)
  )
  
  rate <- rep_len(rate, length.out = n)
  shape <- rep_len(shape, length.out = n)
  theta <- rep_len(theta, length.out = n)
  
  out <- numeric(n)
  i0 <- theta == 0
  i1 <- theta > 0
  
  if (any(i0)) {
    out[i0] <- rate[i0] * shape[i0] * x[i0]^(shape[i0] - 1)
  }
  
  if (any(i1)) {
    out[i1] <- rate[i1] * shape[i1] * x[i1]^(shape[i1] - 1) /
      (1 + theta[i1] * rate[i1] * x[i1]^shape[i1])
  }
  
  out
}

#' @export
#' @rdname GLL
HGLL <- function(x, rate, shape = 1, theta = 0) {
  n <- length(x)
  stopifnot(
    "`shape` must be a positive vector of length 1 or length of `x`." =
      is.numeric(shape) && all(shape > 0) && length(shape) %in% c(1, n),
    "`rate` must be a positive vector of length 1 or length of `x`." =
      is.numeric(rate) && all(rate > 0) && length(rate) %in% c(1, n),
    "`theta` must be a non-negative vector of length 1 or length of `x`." =
      is.numeric(theta) && all(theta >= 0) && length(theta) %in% c(1, n),
    "`x` quantiles must be a non-negative vector." =
      is.numeric(x) && all(x >= 0)
  )
  
  rate <- rep_len(rate, length.out = n)
  shape <- rep_len(shape, length.out = n)
  theta <- rep_len(theta, length.out = n)
  
  out <- numeric(n)
  i0 <- theta == 0
  i1 <- theta > 0
  
  if (any(i0)) {
    out[i0] <- rate[i0] * x[i0]^shape[i0]
  }
  
  if (any(i1)) {
    out[i1] <- log1p(theta[i1] * rate[i1] * x[i1]^shape[i1]) / theta[i1]
  }
  
  out
}


#' @export
#' @rdname GLL
oGLL <- function(x, rate, shape = 1, theta = 0, lower.tail = FALSE) {
  n <- length(x)
  lower_tail <- as.logical(lower.tail[1])

  stopifnot(
    "`shape` must be a positive vector of length 1 or length of `x`." =
      is.numeric(shape) && all(shape > 0) && length(shape) %in% c(1, n),
    "`rate` must be a positive vector of length 1 or length of `x`." =
      is.numeric(rate) && all(rate > 0) && length(rate) %in% c(1, n),
    "`theta` must be a non-negative vector of length 1 or length of `x`." =
      is.numeric(theta) && all(theta >= 0) && length(theta) %in% c(1, n),
    "`x` must be a non-negative numeric vector." =
      is.numeric(x) && all(x >= 0)
  )
  
  rate <- rep_len(rate, length.out = n)
  shape <- rep_len(shape, length.out = n)
  theta <- rep_len(theta, length.out = n)
  
  out <- numeric(n)
  i0 <- theta == 0
  i1 <- theta > 0
  
  if (any(i0)) {
    H <- rate[i0] * x[i0]^shape[i0]
    out[i0] <- 1 / expm1(H)
  }
  
  if (any(i1)) {
    H <- log1p(
      theta[i1] * rate[i1] * x[i1]^shape[i1]
    ) / theta[i1]
    
    # Equivalent to 1 / (exp(H) - 1), but more stable.
    out[i1] <- exp(-H) / (-expm1(-H))
  }
  if(lower_tail){
    out <- 1/out
  }
  out
}