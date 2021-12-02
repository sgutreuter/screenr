#' screenr
#'
#' The \code{screenr} package enables construction of binary test-screening
#' tools.  It is designed to enable those with only a basic familiarity with
#' \code{R} to develop, validate and implement screening tools for diagnostic
#' tests.
#'
#' Consider the situation where a diagnostic test for some condition is
#' relatively expensive, and the condition is rare.  In that case, universal
#' testing would not be efficient in terms of the yield of postive results per
#' test performed.  Now suppose that responses to a set of simple questions may
#' be predictive of the condition.  Package \code{screenr} enables estimation of
#' thresholds for making decisions about when to test in order to screen
#' in/out individuals based on Receiver Operating Characteristics (ROC)
#' estimated from an initial sample.  The choice of a particular screening
#' threshold is left to the user, and should be based on careful consideration
#' of application-specific tradeoffs between sensitivity and specificity.
#' \code{screer} also enables easy construction of screening tools.
#'
#' @details
#' The high-level functions in the \code{screenr} package are:
#' \describe{
#' \item{\code{\link{glmpathScreenr}}}{Test-screening based on GLM path regularization
#' of logistic regression models}
#' \item{\code{\link{logisticScreenr}}}{Test-screening based maximum-likelhood estimation
#' of logistic regression models}
#' \item{\code{\link{easyTool}}}{Easy implementation of test-screening tools}
#' \item{\code{\link{simpleScreenr}}}{Simple un-optimized test-screening}
#' \item{\code{link{rescale_to_Int}}}{Rescale a strictly positive vector of real numbers
#' to integers}
#' \item{\code{\link{sens_spec}}}{Sensitivity and specificity from a 2 x 2 table}
#' }
#'
#' There are \code{plot}, \code{print} \code{summary}, \code{predict},
#' \code{getWhat}, and \code(ntpp) methods for the objects produced by
#' \code{glmpathScreenr}, \code{logisticScreenr}, \code{simpleScreenr}
#' and \code{easyTool}. In addition, there is a \code{coef} method for
#' \code{glmpathScreenr} and \code{logisticScreenr} objects.
#'
#' @note
#' The canonical source for \code{screenr} is
#' \url{https://github.com/sgutreuter/screenr}
#'
#' @name screenr
#' @docType package
#' @author Steve Gutreuter: \email{sgutreuter@gmail.com}
NULL
