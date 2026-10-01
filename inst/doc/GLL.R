## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
R <- function() knitr::include_graphics("Rlogo.png", dpi = 5000)

## ----echo=FALSE, out.width = '33%'--------------------------------------------
knitr::include_graphics("hex-hce.png")

## ----eval = TRUE--------------------------------------------------------------
library(hce)
packageVersion("hce")

## ----eval=FALSE---------------------------------------------------------------
# rweibullGF()
# rGLL()
# pGLL()
# qGLL()
# hGLL()
# HGLL()

## -----------------------------------------------------------------------------
n <- 1000
rate <- 0.5
shape <- 2
theta <- 1.5

y_gll <- rGLL(
  n,
  rate = rate,
  shape = shape,
  theta = theta
)

y_weibull_gf <- rweibullGF(
  n,
  rate = rate,
  shape = shape,
  theta = theta
)

ecdf_gll <- ecdf(y_gll)
ecdf_weibull_gf <- ecdf(y_weibull_gf)

x <- seq(0, 5, 0.01)

plot(
  x,
  ecdf_gll(x),
  type = "s",
  col = "blue",
  lwd = 2,
  lty = 1,
  ylim = c(0, 1),
  ylab = "Empirical Cumulative Distribution Function",
  xlab = "Time",
  main = "Empirical CDFs of GLL and Weibull with Gamma Frailty"
)

lines(
  x,
  ecdf_weibull_gf(x),
  type = "s",
  col = "red",
  lwd = 3,
  lty = 2
)

legend(
  "bottomright",
  legend = c("GLL", "Weibull with gamma frailty"),
  col = c("blue", "red"),
  lwd = c(2, 3),
  lty = c(1, 2)
)

## -----------------------------------------------------------------------------
rate1 <- 0.08
rate0 <- 0.10
shape <- 0.85

hr_at_time_zero <- rate1 / rate0

hazard_ratio <- function(time, theta, alpha = shape) {
  hr_at_time_zero *
    (1 + theta * rate0 * time^alpha) /
    (1 + theta * rate1 * time^alpha)
}

time <- seq(0, 3, 0.001)

plot(
  time,
  hazard_ratio(time, theta = 1),
  type = "l",
  col = "blue",
  lwd = 2,
  ylim = c(0.75, 1.1),
  log = "y",
  ylab = "Hazard Ratio",
  xlab = "Time",
  lty = 1
)

lines(
  time,
  hazard_ratio(time, theta = 5),
  col = "red",
  lwd = 2,
  lty = 2
)

lines(
  time,
  hazard_ratio(time, theta = 10),
  col = "darkgreen",
  lwd = 2,
  lty = 3
)

abline(
  h = c(hr_at_time_zero, 1),
  lty = 4
)

legend(
  "topright",
  legend = expression(theta == 1, theta == 5, theta == 10),
  col = c("blue", "red", "darkgreen"),
  lwd = 2,
  lty = 1:3
)

