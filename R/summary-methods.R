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


## Function summary.easy_tool
##
#' A Summary Method for \code{easy_tool}-Class Objects
#'
#' @param object an \code{easy_tool} object.
#'
#' @param ... optional arguments passed to \code{summary} methods.
#'
#' @return a dataframe containing the summary, including the Df, Deviance,
#' AIC and BIC for each step along the GLM path for which the active set
#' changed.
#'
#' @details This is essentially a wrapper for \code{glmpath::summary.glmpath}
#' provided for \code{lasso_screenr} objects.
#' @examples
#' attach(uniobj1)
#' summary(uniobj1)
#' @import pROC
#' @export
summary.easy_tool <- function(object, ...){
    if(!("easy_tool" %in% class(object)))
        stop("object not easy_tool class")
    print(object$ROC)
}


## Function summary.lasso_screenr
##
#' A Summary Method for \code{lasso_screenr}-Class Objects
#'
#' @param object a lasso_screenr object
#'
#' @param ... optional arguments passed to \code{summary} methods.
#'
#' @return a dataframe containing the summary, including the Df, Deviance,
#' AIC and BIC for each step along the GLM path for which the active set
#' changed.
#'
#' @details This is essentially a wrapper for \code{glmpath::summary.glmpath}
#' provided for \code{lasso_screenr} objects.
#' @examples
#' attach(uniobj1)
#' summary(uniobj1)
#' @export
summary.lasso_screenr <- function(object, ...){
    if(!("lasso_screenr" %in% class(object)))
        stop("object not lasso_screenr class")
    res <- object$Summary
    res
}


## Function summary.logreg_screenr
##
#' A Summary Method for \code{logreg_screenr}-Class Objects
#'
#' @param object an object of class \code{logreg_screenr} produced by function
#' \code{logreg_screenr}.
#'
#' @param diagnostics a logical value; plot model diagnostics if \verb{TRUE}.
#'
#' @param ... optional arguments passed to \code{summary} methods.
#'
#' @return Nothing.  Summaries are printed as a side effect.
#'
#' @examples
#' attach(uniobj2)
#' summary(uniobj2)
#' @export
summary.logreg_screenr <- function(object, ..., diagnostics = FALSE){
    if(!("logreg_screenr" %in% class(object))) stop("object not logreg_screenr class")
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
#' A Summary Method for \code{simple_screenr}-Class Objects
#'
#' @param object an object of class \code{simple_screenr}.
#'
#' @param ... optional arguments passed to \code{summary} methods.
#'
#' @return Nothing.  Thresholds, specificities and sensitivities are printed as
#' a side effect.
#'
#' @examples
#' data(unicorns)
#' toosimple <- simple_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
#'                             data = unicorns)
#' summary(toosimple)
#' @export
summary.simple_screenr <- function(object, ...){
    if(!("simple_screenr" %in% class(object)))
        stop("object not a simple_screenr object")
    cat("Call:\n")
    print(object$Call)
    cat("\nPrevalence (In-sample prevalence of condition):\n")
    print(object$Prevalence)
    cat("\nReceiver Operating Characteristics:\n")
    auc <- round(as.numeric(object$ISroc$auc), digits = 4)
    cat(paste("\nIn-sample area under the ROC curve: ",
              auc, "\n", sep = ""))
}
