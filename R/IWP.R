#' Calculates patient-level individual win proportions
#'
#' @param data a data frame containing subject-level data.
#' @param AVAL a character string specifying the variable containing the ordinal analysis values.
#' @param TRTP a character string specifying the treatment variable.
#' @param ref the reference treatment value.
#' @return the input data frame with rows in their original order and a new column of individual win proportions.
#'  The new column is named using the input `AVAL` value followed by `_`, and two additional columns, `AVAL` and `TRTP`, containing
#' the corresponding analysis values and treatment-group assignments. The
#' `AVAL` and `TRTP` columns facilitate conversion to an `hce` object using
#' [hce::as_hce()].   
#' @export
#' @md
#' @seealso [hce::calcWO()], [hce::calcWO.hce()], [hce::calcWO.formula()].
#' @references Gasparyan SB et al. "Adjusted win ratio with stratification: calculation methods and interpretation." Statistical Methods in Medical Research 30.2 (2021): 580-611. <doi:10.1177/0962280220942558>.
#' @examples
#' # Example 1
#' ## Derive individual win proportions using baseline eGFR
#' dat <- KHCE[, c("TRTPN", "EGFRBL")]
#' dat1 <- IWP(
#'   data = dat,
#'   AVAL = "EGFRBL",
#'   TRTP = "TRTPN",
#'   ref = 2
#' )
#'
#' ## Calculate the win proportion and its standard error
#' WP <- tapply(dat1$EGFRBL_, dat1$TRTPN, mean)
#' VAR <- tapply(
#'   dat1$EGFRBL_,
#'   dat1$TRTPN,
#'   function(x) (length(x) - 1) * var(x) / length(x)
#' )
#' N <- tapply(dat1$EGFRBL_, dat1$TRTPN, length)
#' SE <- sqrt(sum(VAR / N))
#'
#' ## Compare the results with the calcWO() implementation
#' options(digits = 10)
#' c(WP = WP[[1]], SE = SE)
#' calcWO(EGFRBL ~ TRTPN, data = dat, ref = 2)[c("WP", "SE_WP")]
#'
#' ## The output includes standardized AVAL and TRTP columns, allowing it to
#' ## be converted directly to an hce object using as_hce().
#' calcWO(as_hce(dat1))
IWP <- function(data, AVAL, TRTP, ref){
  data <- as.data.frame(data)
  ############# Keep the order of the original data for the output
  data$.IWP_input_order <- seq_len(nrow(data))
  ###############################
  AVAL <- AVAL[1]
  ref <- ref[1]
  TRTP <- TRTP[1]
  
  if (!AVAL %in% base::names(data)) 
    stop(paste0("The variable ", AVAL, " is not in the dataset."))
  if (!TRTP %in% base::names(data)) 
    stop(paste0("The variable ", TRTP, " is not in the dataset."))
  data$AVAL <- data[, AVAL]
  data$TRTP <- data[, TRTP]
  if (length(unique(data$TRTP)) != 2) 
    stop("The dataset should contain two treatment groups.")
  if (!ref %in% unique(data$TRTP)) 
    stop("Choose the reference from the values in TRTP.")
  data$TRTP <- base::ifelse(data$TRTP == ref, "P", "A")
  A <- base::rank(c(data$AVAL[data$TRTP == "A"], data$AVAL[data$TRTP == 
                                                             "P"]), ties.method = "average")
  B <- base::tapply(data$AVAL, data$TRTP, base::rank, ties.method = "average")
  n <- base::tapply(data$AVAL, data$TRTP, base::length)
  n1 <- n[["A"]]
  n0 <- n[["P"]]
  d <- base::data.frame(R1 = A, R2 = base::c(B$A, B$P), TRTP = base::c(base::rep("A", 
                                                                                 n1), base::rep("P", n0)))
  d$R <- d$R1 - d$R2
  d$R0 <- base::ifelse(d$TRTP == "A", d$R/n0, d$R/n1)
  data <- rbind(data[data$TRTP == "A", ], data[data$TRTP == "P", ])
  data[ , paste0(AVAL, "_")] <- d$R0
  ### reorder to keep the original order
  data <- data[order(data$.IWP_input_order), , drop = FALSE]
  data$.IWP_input_order <- NULL
  data
}


