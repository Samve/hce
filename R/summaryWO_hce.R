#' Win odds summary for an hce object
#'
#' @param x an hce object.
#' @param ... additional parameters.
#' @returns a list containing the summary of wins, losses, and ties. It contains the following named objects:
#' * summary a data frame containing number of wins, losses, and ties by treatment group and the overall number of comparisons.
#' * summary_by_GROUP (if `GROUP` variable is specified) a summary data frame by `GROUP`.
#' * WO calculated WO (win odds) and WP (win probability) and their standard errors.
#' @export
#' @md
#' @seealso [hce::calcWO()], [hce::summaryWO()], [hce::summaryWO.data.frame()] methods.
#' @examples
#' dat <- new_hce(HCE4)
#' summaryWO(dat, ref = "P")
summaryWO.hce <- function(x, ...){
  Args <- base::list(...)
  x <- base::as.data.frame(x)
  
  if(!is.null(Args[["ref"]])) ref <- Args[["ref"]]
  else ref <- "P"
  if(! ref %in% c("A", "P")) stop("Choose the reference from the values A or P.")
  
  summaryWO.data.frame(x = x, AVAL = "AVAL", TRTP = "TRTP", ref = ref)
}