#' testscreener
#'
#' The \code{testscreenr} package enables construction of binary test screening
#' tool.  Consider the situation where a diagnostic test for some condition is
#' expensive, and the condition is rare.  In that case, universal testing
#' would not be efficient in terms of the yield of postive results per test
#' performed.  Now suppose that responses to a set of simple questions may be
#' predictive of the condition.  Package \code{testscreenr} enables the
#' thresholds for making decisions about when to test in order to of screen
#' in/out individuals based on Receiver Operating Characteristics (ROC)
#' estimated from an initial sample.  The choice of a particular screening
#' threshold is left to the user, and should be based on careful consideration
#' of application-specific tradeoffs between sensitivity and specificity.
#'
#' @name testscreener
#' @docType package
#' @author Steve Gutreuter: \email{sgutreuter@cdc.gov}
#'
NULL
