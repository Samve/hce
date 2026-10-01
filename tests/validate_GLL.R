library(hce)
########################################
rate <- 0.1
shape <- 0.85
theta <- 0

## The GLL rate parameter corresponds to the following Weibull scale.
scale <- rate^(-1 / shape)

## -------------------------------------------------------------------------
## pGLL() versus stats::pweibull()
## -------------------------------------------------------------------------

q <- c(0, 0.01, 0.1, 1, 3, 10, Inf)

## pGLL: lower-tail probabilities
result <- pGLL(
  q = q,
  rate = rate,
  shape = shape,
  theta = theta,
  lower.tail = TRUE,
  log.p = FALSE
)

expected <- stats::pweibull(
  q = q,
  shape = shape,
  scale = scale,
  lower.tail = TRUE,
  log.p = FALSE
)

stopifnot(isTRUE(all.equal(result, expected)))

## pGLL: upper-tail probabilities
result <- pGLL(
  q = q,
  rate = rate,
  shape = shape,
  theta = theta,
  lower.tail = FALSE,
  log.p = FALSE
)

expected <- stats::pweibull(
  q = q,
  shape = shape,
  scale = scale,
  lower.tail = FALSE,
  log.p = FALSE
)

stopifnot(isTRUE(all.equal(result, expected)))

## pGLL: log lower-tail probabilities
result <- pGLL(
  q = q,
  rate = rate,
  shape = shape,
  theta = theta,
  lower.tail = TRUE,
  log.p = TRUE
)

expected <- stats::pweibull(
  q = q,
  shape = shape,
  scale = scale,
  lower.tail = TRUE,
  log.p = TRUE
)

stopifnot(isTRUE(all.equal(result, expected)))

## pGLL: log upper-tail probabilities
result <- pGLL(
  q = q,
  rate = rate,
  shape = shape,
  theta = theta,
  lower.tail = FALSE,
  log.p = TRUE
)

expected <- stats::pweibull(
  q = q,
  shape = shape,
  scale = scale,
  lower.tail = FALSE,
  log.p = TRUE
)

stopifnot(isTRUE(all.equal(result, expected)))


## -------------------------------------------------------------------------
## qGLL() versus stats::qweibull()
## -------------------------------------------------------------------------

p <- c(0, 0.01, 0.1, 0.5, 0.9, 0.99, 1)

## qGLL: lower-tail probabilities
result <- qGLL(
  p = p,
  rate = rate,
  shape = shape,
  theta = theta,
  lower.tail = TRUE,
  log.p = FALSE
)

expected <- stats::qweibull(
  p = p,
  shape = shape,
  scale = scale,
  lower.tail = TRUE,
  log.p = FALSE
)

stopifnot(isTRUE(all.equal(result, expected)))

## qGLL: upper-tail probabilities
result <- qGLL(
  p = p,
  rate = rate,
  shape = shape,
  theta = theta,
  lower.tail = FALSE,
  log.p = FALSE
)

expected <- stats::qweibull(
  p = p,
  shape = shape,
  scale = scale,
  lower.tail = FALSE,
  log.p = FALSE
)

stopifnot(isTRUE(all.equal(result, expected)))

## qGLL: log lower-tail probabilities
log_p <- log(p)

result <- qGLL(
  p = log_p,
  rate = rate,
  shape = shape,
  theta = theta,
  lower.tail = TRUE,
  log.p = TRUE
)

expected <- stats::qweibull(
  p = log_p,
  shape = shape,
  scale = scale,
  lower.tail = TRUE,
  log.p = TRUE
)

stopifnot(isTRUE(all.equal(result, expected)))

## qGLL: log upper-tail probabilities
result <- qGLL(
  p = log_p,
  rate = rate,
  shape = shape,
  theta = theta,
  lower.tail = FALSE,
  log.p = TRUE
)

expected <- stats::qweibull(
  p = log_p,
  shape = shape,
  scale = scale,
  lower.tail = FALSE,
  log.p = TRUE
)

stopifnot(isTRUE(all.equal(result, expected)))


## -------------------------------------------------------------------------
## Vector-valued parameter validation
## -------------------------------------------------------------------------

q <- c(0, 0.1, 1, 10)

rate <- c(0.1, 0.2, 0.5, 1)
shape <- c(0.5, 0.85, 1.5, 2)
theta <- rep(0, length(q))

## Calculate the equivalent Weibull scale for each rate-shape pair.
scale <- rate^(-1 / shape)

## pGLL with vector-valued rate and shape
for (lower_tail in c(TRUE, FALSE)) {
  for (log_p in c(TRUE, FALSE)) {
    result <- pGLL(
      q = q,
      rate = rate,
      shape = shape,
      theta = theta,
      lower.tail = lower_tail,
      log.p = log_p
    )
    
    expected <- stats::pweibull(
      q = q,
      shape = shape,
      scale = scale,
      lower.tail = lower_tail,
      log.p = log_p
    )
    
    stopifnot(isTRUE(all.equal(result, expected)))
  }
}

## qGLL with vector-valued rate and shape
p <- c(0.01, 0.1, 0.5, 0.99)

for (lower_tail in c(TRUE, FALSE)) {
  for (log_p in c(TRUE, FALSE)) {
    p_input <- if (log_p) log(p) else p
    
    result <- qGLL(
      p = p_input,
      rate = rate,
      shape = shape,
      theta = theta,
      lower.tail = lower_tail,
      log.p = log_p
    )
    
    expected <- stats::qweibull(
      p = p_input,
      shape = shape,
      scale = scale,
      lower.tail = lower_tail,
      log.p = log_p
    )
    
    stopifnot(isTRUE(all.equal(result, expected)))
  }
}

cat("All theta = 0 pGLL() and qGLL() validation tests passed.\n")