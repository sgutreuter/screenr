#' screenr
#'
#' The \code{screenr} package enables construction of binary test-screening
#' tools.  Consider the situation where a diagnostic test for some condition is
#' expensive, and the condition is rare.  In that case, universal testing
#' would not be efficient in terms of the yield of postive results per test
#' performed.  Now suppose that responses to a set of simple questions may be
#' predictive of the condition.  Package \code{screenr} enables estimation of
#' thresholds for making decisions about when to test in order to screen
#' in/out individuals based on Receiver Operating Characteristics (ROC)
#' estimated from an initial sample.  The choice of a particular screening
#' threshold is left to the user, and should be based on careful consideration
#' of application-specific tradeoffs between sensitivity and specificity.
#'
#' @name screener
#' @docType package
#' @author Steve Gutreuter: \email{sgutreuter@@gmail.gov}
#' @seealso
#' The high-level functions in the \code{screenr} package are:
#' \describe{
#' \item{\code{binomialScreening}}{Test-Screening Tool Based on Binomial Regression}
#' \item{\code{inverseLink}}{Inverse of Binomial Regression Link Functions}
#' \item{\code{nSensSpec}}{Sample Size for Joint Testing of Sensitivity and Specificity}
#' \item{\code{sens_spec}}{Sensitivity and Specificity from a 2 x 2 Table}
#' \item{\code{simpleScreening}}{Simple Un-optimized Test-Screening Tool}
#' \item{\code{testCounts}}{Expected Number of Tests Required per Positive Test Result}
#'}
#'
#' In addition, there are \code{plot}, \code{print} and \code{summary} methods
#' for the objects produced by \code{binomialScreening} and \code{simpleScreening}.
#'
NULL
