#' UIV Testing Data on Unicorns
#'
#' A preliminary study was conducted in which a random sample of 5,000 properly
#' consented unicorns were recruited from 20 clinics.  Each unicorn was asked
#' five questions about their behavior and health.  Unicorns responded by
#' stomping a hoof once to indicate "no", and twice to indicate "yes".  A
#' sample of venous blood was drawn from each, and was subsequently tested
#' for the presence of antibodies to Unicorn Immunodeficiency Virus (UIV) using
#' a standard assay algorithm.
#'
#' @format A data frame with six columns:
#' \describe{
#' \item(\code{clinic}){Clinic code}
#' \item{\code{Q1}}{Response to screening question 1 (0 = "no", 1 = "yes")}
#' \item{\code{Q2}}{Response to screening question 2 (0 = "no", 1 = "yes")}
#' \item{\code{Q3}}{Response to screening question 3 (0 = "no", 1 = "yes")}
#' \item{\code{Q4}}{Response to screening question 4 (0 = "no", 1 = "yes")}
#' \item{\code{Q5}}{Response to screening question 5 (0 = "no", 1 = "yes")}
#' \item{\code{testresult}}{UIV status, where 0 and 1 denote negative and positive test results, repectively.}
#' }
#' @examples
#' \dontrun{
#' unicorns
#' }
"unicorns"
