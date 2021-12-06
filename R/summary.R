################################################################################
##     R Script: summary.R
##
##      Package: screenr
##
##  Description: S3 methods for summary of screenr objects
##
##       Author: Steve Gutreuter
##               sgutreuter@gmail.com
################################################################################


## Function summary.glmpathScreenr
##
#' \code{summary.glmpathScreenr} returns a summary of the GLM path regularizer
#'
#' @param object a glmpathScreenr object
#'
#' @return a dataframe containing the summary, including the Df, Deviance,
#' AIC and BIC for each step along the GLM path for which the active set
#' changed.
#'
#' @details This is essentially a wrapper for \code{glmpath::summary.glmpath}
#' provided for \code{glmpathScreenr} objects.
#' @examples
#' attach(uniobj1)
#' summary(uniobj1)
#' @export
summary.glmpathScreenr <- function(object){
    if(!("glmpathScreenr" %in% class(object)))
        stop("object not glmpathScreenr class")
    res <- object$Summary
    res
}


## Function summary.logisticScreenr
##
#' \code{summary.logisticScreenr} is an S3 print method for \code{logisticScreenr}
#' objects.
#'
#' @param object an object of class \code{logisticScreenr} produced by function
#' \code{logisticScreenr}.
#'
#' @param diagnostics a logical value; plot model diagnostics if \verb{TRUE}.
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @return Nothing.  Summaries are printed as a side effect.
#'
#' @examples
#' load(uniobj2)
#' class(uniobj2)
#' summary(uniobj2)
#' @export
summary.logisticScreenr <- function(object, diagnostics = FALSE, ...){
    if(!("logisticScreenr" %in% class(object))) stop("object not logisticScreenr class")
    cat("Call:\n")
    print(object$Call)
    cat("\n\nLogistic regression model summary:\n")
    print(summary(object$ModelFit))
    if(diagnostics) plot(object$ModelFit)
    cat("\nPrevalence (In-sample prevalence of condition):\n")
    print(object$Prevalence)
    cat("\nReceiver Operating Characteristics:\n")
    is.auc <- round(as.numeric(object$ISroc$auc), digits = 4)
    cat("\nIn-sample (overly optimistic) area under the curve: ",
        is.auc, "\n", sep = "")
    cv.auc <- round(as.numeric(object$CVroc$auc), digits = 4)
    cat("Out-of-sample area under the curve: ",
        cv.auc, "\n", sep = "")
}


## Function summary.simpleSreenr
##
#' \code{summary.simpleScreenr} is a summary method for \code{simpleScreenr} objects
#'
#' @param object an object of class \code{simpleScreenr}.
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @return Nothing.  Thresholds, specificities and sensitivities are printed as
#' a side effect.
#'
#' @examples
#' data(unicorns)
#' toosimple <- simpleScreenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6,
#'                           data = unicorns)
#' summary(toosimple)
#' @export
summary.simpleScreenr <- function(object, ...){
    if(!("simpleScreenr" %in% class(object)))
        stop("object not a simpleScreenr object")
    cat("Call:\n")
    print(object$Call)
    cat("\nPrevalence (In-sample prevalence of condition):\n")
    print(object$Prevalence)
    cat("\nReceiver Operating Characteristics:\n")
    auc <- round(as.numeric(object$ISroc$auc), digits = 4)
    cat(paste("\nIn-sample area under the ROC curve: ",
              auc, "\n", sep = ""))
}
