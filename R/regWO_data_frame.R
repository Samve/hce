#' Win odds regression using a data frame
#'
#' @param x a data frame containing subject-level data.
#' @param AVAL variable in the data with ordinal analysis values.
#' @param TRTP the treatment variable in the data.
#' @param COVAR a numeric covariate.
#' @param ref the reference treatment group.
#' @param alpha significance level. The default is 0.05.
#' @param WOnull the null hypothesis. The default is 1.
#' @param ... additional parameters.
#' @returns a data frame containing the win odds and its confidence interval. The data frame has an attribute called "covar_info" giving summary statistics for the covariate used for the calculations. The data frame itselfs contains the following columns:
#' * WO_beta adjusted win odds.
#' * LCL lower confidence limit for adjusted WO.
#' * UCL upper confidence limit for adjusted WO.
#' * SE standard error of the adjusted win odds.
#' * WOnull win odds of the null hypothesis (specified in the `WOnull` argument).
#' * alpha two-sided significance level for calculating the confidence interval (specified in the `alpha` argument).
#' * Pvalue p-value associated with testing the null hypothesis.
#' * beta adjusted win probability.
#' * SE_beta standard error for the adjusted win probability.
#' * WP (non-adjusted) win probability.
#' * SE_WP standard error of the non-adjusted win probability.
#' * WO non-adjusted win odds.
#' @export
#' @md
#' @seealso [hce::regWO()].
#' @references Gasparyan, Samvel B., et al. "Adjusted win ratio with stratification: calculation methods and interpretation." Statistical Methods in Medical Research 30.2 (2021): 580-611. <doi:10.1177/0962280220942558>
#' @examples
#' # An baseline covariate that is highly correlated with the outcome
#' set.seed(2023)
#' dat <- COVID19
#' n <- nrow(dat)
#' dat$Severity <- ifelse(dat$GROUP > 4, rexp(n, 1), rexp(n, 10))
#' res <- regWO(x = dat, AVAL = "GROUP", TRTP = "TRTP", COVAR = "Severity", ref = "Placebo")
#' res
#' attr(res, "covar_info")
regWO.data.frame <- function(x, AVAL, TRTP, COVAR, ref, alpha = 0.05, WOnull = 1, ...){
  data <- as.data.frame(x)
  alpha <- alpha[1]
  WOnull <- WOnull[1]
  WPnull <- WOnull/(WOnull + 1)
  Ca <- stats::qnorm(1 - alpha/2)
  
  data$AVAL <- data[, base::names(data) == AVAL]
  data$TRTP <- data[, base::names(data) == TRTP]
  data$COVAR <- data[, base::names(data) == COVAR]
  if(length(unique(data$TRTP)) != 2) stop("The dataset should contain two treatment groups.")
  if(!ref %in% unique(data$TRTP)) stop("Choose the reference from the values in TRTP.")
  data$TRTP <- base::ifelse(data$TRTP == ref, "P", "A")
  if(!is.numeric(data$COVAR)) stop("COVAR should be numeric.")  
 
  A <- base::rank(c(data$AVAL[data$TRTP == "A"], data$AVAL[data$TRTP == "P"]), ties.method = "average")
  B <- base::tapply(data$AVAL, data$TRTP, base::rank, ties.method = "average")
  n <- base::tapply(data$AVAL, data$TRTP, base::length)
  n1 <- n[["A"]]
  n0 <- n[["P"]]
  N <- n0 + n1
  d <- base::data.frame(R1 = A, R2 = base::c(B$A, B$P), 
                        TRTP = base::c(base::rep("A", n1), base::rep("P", n0)), 
                        COVAR = data$COVAR)
  d$R <- d$R1 - d$R2
  d$R0 <- base::ifelse(d$TRTP == "A", d$R/n0, d$R/n1)
  
  WP0 <- base::tapply(d$R0, d$TRTP, base::mean)
  WP <- WP0[["A"]]
  VAR <- base::tapply(d$R0, d$TRTP, function(x) (base::length(x)-1)*stats::var(x)/base::length(x))
  SE_WP <- base::sqrt(base::sum(VAR/n))
  
  
  MEAN_COVAR <- base::tapply(d$COVAR, d$TRTP, base::mean)
  VAR_COVAR <- base::tapply(d$COVAR, d$TRTP, function(x) (base::length(x)-1)*stats::var(x)/base::length(x))
  C0 <- sapply(split(d[, c("COVAR", "R0")], d$TRTP), 
               function(d) (length(d$COVAR) - 1)*stats::cov(d$COVAR, d$R0)/length(d$COVAR))
  C1 <- sum(C0/c(n1, n0))
  
  beta <- WP - C1*(MEAN_COVAR["A"] - MEAN_COVAR["P"])/sum(VAR_COVAR/c(n1, n0))
  SE_beta2 <-  SE_WP^2 - C1^2/sum(VAR_COVAR/c(n1, n0))
  SE_beta <- sqrt(SE_beta2)
  
  WO_beta <- beta/(1 - beta)
  SE <- SE_beta/(beta*(1 - beta))
  LCL <- WO_beta*base::exp(- Ca*SE)
  UCL <- WO_beta*base::exp(Ca*SE)
  threshold <- base::abs(beta - WPnull)/SE_beta
  P <- 2*(1 - stats::pnorm(threshold))
  
  out <- base::data.frame(WO_beta = WO_beta, LCL = LCL, UCL = UCL, 
                          SE = SE, WOnull = WOnull, 
                          alpha = alpha, Pvalue = P, 
                          beta = beta, SE_beta = SE_beta, 
                          WP = WP, SE_WP = SE_WP, WO = WP/(1 - WP))
  
  covar_info <- data.frame(MEAN_DIFF = MEAN_COVAR["A"] - MEAN_COVAR["P"], VAR_SUM = sum(VAR_COVAR/c(n1, n0)), COV_SUM = C1)
  attr(out, "covar_info") <- covar_info
  
  return(out)
}
