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
##
#################################################################################


#' A Simple Un-optimized HIV Testing Screener
#'
#' \code{simpleScreening} Development of a simple screening tool.
#'
#' Development of a simple testing screening tool based on the sum of ordinal
#' numeric responses to questions which are predictive of the testing results.
#' \code{binomialScreening} will generally outperform this approach, and
#' \code{simpleScreening} is provided mainly for comparison.
#'
#' S3 \code{plot}, \code{print} and \code{summary} methods are available.
#'
#' @param formula A formula for the screener.
#' @param data A data frame containing the testing outcome and predictive
#' covariates to be used for testing screening.  The testing outcome must
#' be binary (0,1) indicating negative and positive test results, respectively.
#' The covariates are typically binary (0 = no, 1 = yes) responses to
#' questions which may be predictive of the test result.
#'
#' @return A list having S3 class "simplescreenr" containing the call, the
#' Receiver Operating Characteristics and the scores.
#'
#' @author Steve Gutreuter, \email{sgutreuter@@gmail.com}
#'
#' @references Bandason et al 2016. Validation of a screening tool to
#' identify older children living with HIV in primary care facilities in high
#' HIV prevalence settings. AIDS 30(5):779-785
#' \url{http://dx.doi.org/10.1097/QAD.0000000000000959
#'
#' @examples
#' data(unicorns)
#' simpletool <- simpleScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 ,
#'                        data = unicorns)
#' summary(simpletool)
#' plot(simpletool)
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
    scores <- cbind(dat, Bandason.score = score)
    result <- list(Call = call, roc = roc, scores = scores)
    class(result) <- "simplescreenr"
    invisible(result)
}


#' An S3 Print Method
#'
#' \code(print.simplescreenr) An S3 printing method for objects created by
#' function \code(simpleScreening()).
#'
#' @param x A simplescreenr class object.
print.simplescreenr <- function(x){
    if(!class(x) == "simplescreenr") stop("x not a simplescreenr object")
    x$roc
}


#' An S3 Plot Method
#'
#' \code(plot.simplescreenr) An S3 plotting method for objects created by
#' function \code(simpleScreening()).
#'
#' @param x A simplescreenr class object.
#' @export
plot.simplescreenr <- function(x){
    if(!("simplescreenr" %in% class(x))) stop("x must be a simplescreenr object")
    plot.ROC(x$roc)
}


#' An S3 Summary Method
#'
#' \code{summary.simplescreenr} An S3 summary method for objects created by
#' function \code{simpleScreening()}.
#'
#' @param x A simplescreenr class object.
#' @export
summary.simplescreenr <- function(x){
    if(!class(x) == "simplescreenr") stop("x not a Bandason object")
    cat("A list object with three components:\n")
    cat("\nCall:\n")
    print(x$Call)
    cat("\nroc (In-sample Receiver Operating Characteristic):\n")
    print(x$roc)
}

################################   END of FILE   ################################
