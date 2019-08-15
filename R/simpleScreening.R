#################################################################################
##       R PROGRAM: simpleScreening.R
##
##         PROJECT: R fuctions for HIV screening tool development
##
##      WRITTEN BY: Steve Gutreuter, CDC/CGH/DGHT/Statistics, Estimation and
##                               Modeling Team
##                  E-mail:  sgutreuter@cdc.gov
##
##      DISCLAIMER: Although this code has been used by the Centers for Disease
##                  Control & Prevention (CDC), no warranty, expressed or
##                  implied, is made by the CDC or the U.S. Government as to the
##                  accuracy and functioning of the code and related program
##                  material nor shall the fact of distribution constitute any
##                  such warranty, and no responsibility is assumed by the CDC in
##                  connection therewith.
#################################################################################

#' Simple Un-optimized Test-Screening Tool
#'
#' Compute the in-sample performances for development of a very simple test
#' screening tool. The results provide information from which to choose a
#' counting-based threshold score above which a diagnostic test would be
#' performed. \code{binomialScreening} will almost certainly outperform this
#' approach, and \code{simpleScreening} is intended for use in those situations
#' where post-estimation technical capacity is limited to counting responses to
#' questions.
#'
#' @param formula an object of class \code{\link{formula}} defining the testing
#' outcome and predictor covariates.
#' @param data the "training" sample; a data frame containing the testing outcome
#' and predictive covariates to be used for testing screening.  The testing
#' outcome must be binary (0,1) indicating negative and positive test results,
#' respectively, or logical (TRUE/FALSE).  The covariates are typically binary
#' (0 = no, 1 = yes) responses to questions, but the responses may also be
#' ordinal numeric values.
#'
#' @return An object of class "simplescreenr" containing the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{Prevalence}}{Prevalence of the test condition in the training sample.}
#' \item{\code{InSamplePerf}}{A data frame containing in-sample (overly optimistic)
#' sensitivities and specificities.}
#' \item{\code{Scores}}{The training sample, including the scores.}
#' }
#'
#' @references Bandason et al 2016. Validation of a screening tool to
#' identify older children living with HIV in primary care facilities in high
#' HIV prevalence settings. AIDS 30(5):779-785
#' \url{http://dx.doi.org/10.1097/QAD.0000000000000959}
#'
#' @examples
#' data(unicorns)
#' simple <- simpleScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 ,
#'                        data = unicorns)
#' summary(simple)
#' plotROC(simple)
#' testCounts(simple)
#'
#' @export
simpleScreening <- function(formula, data){
    mf <- match.call(expand.dots = FALSE)
    call <- mf
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    dat <- data[complete.cases(data[, names(mf), ]), ]
    y <- model.response(mf, "numeric")
    if(!all(y %in% c(0, 1))) stop("Response variable must be binary (0, 1)")
    prev <- mean(y, na.rm = TRUE)
    preds <- mf[, -1]
    npreds <- dim(preds)[2]
    score <- apply(preds, 1, sum)
    critvals <- 1L:npreds
    tst <- outer(score, critvals, `>=`)
    roc <- data.frame(score = NULL, sensitivity = NULL, specificity = NULL)
    for(i in seq_along(critvals)){
        tabl <- table(factor(tst[, i] * 1 , levels = c("0", "1")),
                      factor(y, levels = c("0", "1")))
        ss <- sens_spec(tabl)
        roc <- rbind(roc, data.frame(score = i,
                                     sensitivity = ss$sensitivity,
                                     specificity = ss$specificity))
        class(roc) <- c("ROC", "data.frame")
    }
    scores <- cbind(dat, score = score)
    result <- list(Call = call,
                   Prevalence = prev,
                   InSamplePerf = roc,
                   Scores = scores)
    class(result) <- "simplescreenr"
    invisible(result)
}


#' An S3 Print Method for \code{simplescreenr} Objects
#'
#' @param x a \code{simplescreenr} class object.
#' @param ... further arguments passed to or from other methods.
#' @param quote logical, indicating whether or not strings should be printed
#' with surrounding quotes.
#' @export
print.simplescreenr <- function(x, quote = FALSE, ...){
    if(!("simplescreenr" %in% class(x))) stop("x not a simplescreenr object")
    x$InSamplePerf
}


#' An S3 Summary Method for \code{simplescreenr} Objects
#'
#' @param object A \code{simplescreenr} class object.
#' @param ... further arguments passed to or from other methods.
#' @export
summary.simplescreenr <- function(object, ...){
    if(!("simplescreenr" %in% class(object))) stop("object not a simplescreenr object")
    cat("\nCall:\n")
    print(object$Call)
    cat("\nPrevalence (In-sample prevalence of condition)\n")
    print(object$Prevalence)
    cat("\nInSampPerf (In-sample Receiver Operating Characteristics):\n")
    print(object$InSamplePerf)
}

################################   END of FILE   ################################
