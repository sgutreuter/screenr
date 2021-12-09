#################################################################################
##       R PROGRAM: simpleScreening.R
##
##         PROJECT: Screening based on un-optimized response counts
##
##      WRITTEN BY: Steve Gutreuter
##                  E-mail:  sgutreuter@gmail.com
##
#################################################################################

## Function simpleScreenr
##
#' \code{simpleScreenr} is used to produce overly simple test-screening tools.
#'
#' @param formula an object of class \code{\link[stats]{formula}} defining the
#' testing outcome and predictor covariates.
#'
#' @param data the "training" sample; a data frame containing the testing outcome
#' and predictive covariates to be used for testing screening.  The testing
#' outcome must be binary (0,1) indicating negative and positive test results,
#' respectively, or logical (\verb{TRUE}/\verb{FALSE}), and the screening scores
#' are the row-wise sums of the values of those covariates.  The covariates are
#' typically binary (0 = no, 1 = yes) responses to questions, but the responses
#' may also be ordinal numeric values.
#'
#' @return An object of class \code{simpleScreenr} containing the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{Prevalence}}{Prevalence of the test condition in the training sample.}
#' \item{\code{ISroc}}{An object of class \code{\link[pROC]{roc}} containing
#' the "in-sample" (overly-optimistic) receiver operating characteristics,
#' and additional functions for use with this object are available in the
#' \code{pROC} package.}
#' \item{\code{Scores}}{The training sample, including the scores.}
#' }
#'
#' @details
#' \code{simpleScreenr} computes the in-sample (\emph{overly optimistic})
#' performances for development of a very simple test screening tool implementing
#' the method of Bandason et al. (2016).  \code{simpleScreener} is not optimized
#' and is intended only for comparision with \code{glmpathScreenr} or
#' \code{logisticScreenr}, either of which will almost certainly  out-perform
#' \code{simpleScreening}.
#'
#' @references
#' Bandason T, McHugh G, Dauya E, Mungofa S, Munyati SM, Weiss HA, Mujuru H,
#' Kranzer K, Ferrand RA. Validation of a screening tool to
#' identify older children living with HIV in primary care facilities in high
#' HIV prevalence settings. AIDS. 2016;30(5):779-785
#' \url{http://dx.doi.org/10.1097/QAD.0000000000000959}
#'
#' Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez J-C,
#' Müller M. \code{pROC}: An open-source package for \code{R} and S+ to
#' analyze and compare ROC curves. BMC Bioinformatics. 2011;12(77):1-8.
#' \url{http://doi.org/10.1186/1471-2105-12-77}
#'
#'
#' @examples
#' data(unicorns)
#' toosimple <- simpleScreenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
#'                            data = unicorns)
#' summary(toosimple)
#'
#' @seealso \code{\link{glmpathScreenr}}, \code{\link{logisticScreenr}}
#' @import pROC
#' @importFrom stats model.response
#' @export
simpleScreenr <- function(formula, data){
    warning("WARNING! WARNING! WARNING! simpleScreenr is suboptimal and is provided only for comparison with other methods." )
    mf <- match.call(expand.dots = FALSE)
    call <- mf
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    dat <- data[stats::complete.cases(data[, names(mf), ]), ]
    y <- stats::model.response(mf, "numeric")
    if(!all(y %in% c(0, 1))) stop("Response variable must be binary (0, 1)")
    prev <- mean(y, na.rm = TRUE)
    preds <- mf[, -1]
    npreds <- dim(preds)[2]
    score <- apply(preds, 1, sum)
    is.roc <- pROC::roc(y, score, auc = TRUE, direction = "<")
    scores <- cbind(dat, score = score)
    result <- list(Call = call,
                   Prevalence = prev,
                   ISroc = is.roc,
                   Scores = scores)
    class(result) <- "simpleScreenr"
    invisible(result)
}
