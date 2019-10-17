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
#' Compute the in-sample (\emph{overly optimistic}) performances for development of a
#' very simple test screening tool. The results provide information from which to
#' choose a counting-based threshold score above which a diagnostic test would be
#' performed. However, it is unlikely that the target sensitivity and specificity
#' would be acheived from new data.  \code{simpleScreening} is intended for use
#' \emph{only} in those situations where the technical capacity for implementation
#' of screening is limited to counting responses to questions.  The methods
#' implemented by \code{binomialScreening} and \code{mebinomScreening} will almost
#' certaily out-perform \code{simpleScreening}.
#'
#' @param formula an object of class \code{\link[stats]{formula}} defining the testing
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
#' plot(simple)
#' \dontrun{testCounts(simple)}
#'
#' @seealso \code{\link{binomialScreening}}, \code{\link{mebinomScreening}}
#' @import pROC
#' @importFrom stats model.response complete.cases
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
    is.roc <- pROC::roc(y, score, auc = TRUE, ci = TRUE, of = "sp",
                        se = seq(0, 1, 0.025), ci.type = "shape")
    scores <- cbind(dat, score = score)
    result <- list(Call = call,
                   Prevalence = prev,
                   ISroc = is.roc,
                   Scores = scores)
    class(result) <- "simplescreenr"
    invisible(result)
}


#' Print Summaries of \code{simplescreenr} Objects
#'
#' @param object A \code{simplescreenr} class object.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Nothing.  Thresholds, specificities and sensitivities are printed as
#' a side effect.
#'
#' @seealso \code{\link{getROC}}
#' @export
summary.simplescreenr <- function(object, ...){
    if(!("simplescreenr" %in% class(object)))
        stop("object not a simplescreenr object")
    cat("Call:\n")
    print(object$Call)
    cat("\nPrevalence (In-sample prevalence of condition):\n")
    print(object$Prevalence)
    cat("\nReceiver Operating Characteristics:\n")
    auc <- round(as.numeric(object$ISroc$auc), digits = 4)
    cat(paste("\nIn-sample (overly optimistic) area under the curve: ",
              auc, "\n", sep = ""))
}

#' Plot ROC Curves of \code{simplescreenr} Objects
#'
#' Plot cross-validated (out-of-sample) ROC curve with pointwise 95% confidence
#' intevals on specificity (gray shaded region), along with the overly optimistic
#' in-sample ROC curve.
#'
#' @param x A object of class "simplescreenr".
#' @param ... Additional arguments passed to \code{\link{plot.roc}} and friends.
#'
#' @return Nothing.  This function produces a plot as a side effect
#'
#' @references
#' Fawcett T. An introduction to ROC analysis. Pattern Recognition Letters. 2006.
#' 27(8):861-874.
#' \url{https://doi.org/10.1016/j.patrec.2005.10.010}
#'
#' Linden A. Measuring diagnostic and predictive accuracy in disease
#' management: an introduction to receiver operating characteristic (ROC) analysis.
#' Journal of Evaluation in Clinical Practice. 2006; 12(2):132-139.
#' \url{https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1365-2753.2005.00598.x}
#'
#' Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez J-C, Muller M.
#' pROC: an open-source package for R and S+ to analyze and compare ROC curves.
#' BMC Bioinformatics 2011; 12:77. \url{https://www.biomedcentral.com/1471-2105/12/77}
#' @importFrom graphics plot
#' @export
plot.simplescreenr <- function(x, ...){
    if(!class(x) == "simplescreenr") stop("x is not a simplescreenr object")
    rocobj1 <- plot(x$ISroc, print.auc = TRUE, ci = TRUE, of = "sp",
                    se = seq(0, 1, 0.025), ci.type = "shape", ...)
}


#' Print Receiver Operating Characteristics for \code{simplescreenr} Objects
#'
#' @param x a \code{simplescreenr} class object.
#' @param ... further arguments passed to or from other methods.
#' @param quote logical, indicating whether or not strings should be printed
#' with surrounding quotes.
#'
#' @return Nothing. Thresholds, specificities and sensitivities are printed as a
#' side effect.
#'
#' @references
#' Fawcett T. An introduction to ROC analysis. Pattern Recognition Letters. 2006.
#' 27(8):861-874.
#' \url{https://doi.org/10.1016/j.patrec.2005.10.010}
#'
#' Linden A. Measuring diagnostic and predictive accuracy in disease
#' management: an introduction to receiver operating characteristic (ROC) analysis.
#' Journal of Evaluation in Clinical Practice. 2006; 12(2):132-139.
#' \url{https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1365-2753.2005.00598.x}
#'
#' Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez J-C, Muller M.
#' pROC: an open-source package for R and S+ to analyze and compare ROC curves.
#' BMC Bioinformatics 2011; 12:77. \url{https://www.biomedcentral.com/1471-2105/12/77}
#'
#' @seealso \code{\link{getROC}}
#' @export
print.simplescreenr <- function(x, quote = FALSE, ...){
    if(!("simplescreenr" %in% class(x))) stop("x not a simplescreenr object")
    cat("\nIn-sample (overly optimistic) sensitivity and specificity:\n")
    df_ <- data.frame(score = ceiling(x$ISroc$thresholds),
                      sensitivity = x$ISroc$sensitivities,
                      specificity = x$ISroc$specificities)
    print(df_)
}

################################   END of FILE   ################################
