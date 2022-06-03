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
#' An S3 Summary Method for \code{screenr} Objects
#'
#' @description \code{summary.easy_tool} provides a summary method for
#' \code{easy-tool}-class objects.
#'
#' @param object an \code{easy_tool} object.
#'
#' @param ... optional arguments passed to \code{summary} methods.
#'
#' @return Nothing.  A summary is printed as a side effect.
#'
#' @examples
#' attach(uniobj1)
#' summary(uniobj1)
#' @import pROC
#' @export
summary.easy_tool <- function(object, ...){
    if(!("easy_tool" %in% class(object)))
        stop("object not easy_tool class")
    cat("\nCall:\n")
    print(object$Call)
    cat("\nRescaled question weights:\n")
    print(object$QuestionWeights )
    cat("\nArea under the ROC curve:\n")
    print(object$ROC)
}


## Function summary.lasso_screenr
##
#' An S3 Summary Method for \code{screenr} Objects
#'
#' @description \code{summary.lasso_screenr} provides a summary method for
#' \code{lasso_screenr}-class objects.
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
#' An S3 Summary Method for \code{screenr} Objects
#'
#' @description \code{summary.logreg_screenr} provides a summary method for
#' \code{logreg_screenr}-class objects.
#'
#' @param object an object of class \code{logreg_screenr} produced by function
#' \code{logreg_screenr}.
#'
#' @param diagnostics a logical value; plot model diagnostics if \verb{TRUE}.
#'
#' @param ... optional arguments passed to \code{summary} methods.
#'
#' @return Nothing. A summary is printed as a side effect.
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


## Function summary.gee_screenr
##
#' An S3 Summary Method for \code{screenr} Objects
#'
#' @description \code{summary.gee_screenr} provides a summary method for
#' \code{gee_screenr}-class objects.
#'
#' @param object an object of class \code{gee_screenr} produced by function
#' \code{gee_screenr}.
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
summary.gee_screenr <- function(object, ..., diagnostics = FALSE){
    if(!("gee_screenr" %in% class(object))) stop("object not gee_screenr class")
    cat("Call:\n")
    print(object$Call)
    cat("\n\nLogistic regression model summary:\n")
    print((summary(object$ModelFit))$coefficients)
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
#' An S3 Summary Method for \code{screenr} Objects
#'
#' @description \code{summary.simple_screenr} provides a summary method for
#' \code{simple_screenr}-class objects.
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
