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
#' @details
#' The high-level functions in the \code{screenr} package are:
#' \describe{
#' \item{\code{\link{binomialScreening}}}{Test-screening based on binomial regression}
#' \item{\code{\link{getCoefficients}}}{Extract the estimated coefficients from \code{binomscreenr} objects}
#' \item{\code{\link{getOddsRatios}}}{Extract a dataframe containing odds ratios and their profile-likelihood confidence intervals}
#' \item{\code{\link{getROC}}}{Extract Receiver Operating Characteristics from objects}
#' \item{\code{\link{inverseLink}}}{Inverse of binomial regression link functions}
#' \item{\code{\link{mebinomScreening}}}{Test-screening based on mixed-effects binomial regression
#' produced by \code{binomialScreening}, \code{mebinomScreening} and \code{simpleScreening}}
#' \item{\code{\link{sens_spec}}}{Sensitivity and specificity from a 2 x 2 table}
#' \item{\code{\link{simpleScreening}}}{Simple un-optimized test-screening}
#' \item{\code{\link{testCounts}}}{Expected number of tests required per positive test result}
#' }
#'
#' In addition, there are \code{plot}, \code{print} and \code{summary} methods
#' for the objects produced by \code{binomialScreening}, \code{mebinomScreening} and
#' \code{simpleScreening}.
#'
#' @name screenr
#' @docType package
#' @author Steve Gutreuter: \email{sgutreuter@@gmail.com}
NULL
