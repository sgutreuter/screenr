################################################################################
##     R Script: print.R
##
##      Package: screenr
##
##  Description: S3 methods for printing screenr objects
##
##       Author: Steve Gutreuter
##               sgutreuter@gmail.com
################################################################################


## Function print.easy_tool
##
#' An S3 Print Method for \code{screenr} Objects
#'
#' @description \code{print.easy_tool} is a print method.
#'
#' @param x an object of class \code{easy_tool}.
#'
#' @param ... optional arguments to \code{print} methods.
#'
#' @examples
#' attach(uniobj1)
#' print(uniobj1)
#' @export
print.easy_tool <- function(x, ...){
    cat("Function call:\n")
    print(x$Call)
    cat("\nType:\n")
    print(x$Type )
    cat("\nQuestion weights:\n")
    print(x$QuestionWeights )
}


## Function print.lasso_screenr
##
#' An S3 Print Method for \code{screenr} Objects
#'
#' @description \code{print.lasso_screenr} is a print method for
#' \code{lasso_screenr}-class objects.
#'
#' @param x an object of class \code{lasso_screenr}
#'
#' @param ... optional arguments to \code{print} methods.
#'
#' @examples
#' attach(uniobj1)
#' print(uniobj1)
#' @export
print.lasso_screenr <- function(x, ...){
    cat("Function call:\n")
    print(x$Call)
    cat("\nglmpath x:\n")
    glmpath::print.glmpath(x$glmpathObj)
    cat("\nTest positivity:", x$Prevalence, "\n")
}


## Function print.logreg_screenr
##
#' An S3 Print Method for \code{screenr} Objects
#'
#' @description \code{print.logreg_screenr} is a print method for
#' \code{logreg_screenr}-class objects.
#'
#' @param x an object of class \code{logreg_screenr}.
#'
#' @param quote logical indicator for whether or not strings should be printed.
#'
#' @param ... optional arguments to \code{print} methods.
#'
#' @return
#' Nothing. Thresholds, specificities and sensitivities are printed as a
#' side effect.
#'
#' @examples
#' attach(uniobj2)
#' print(uniobj2)
#'
#' @export
print.logreg_screenr <- function(x, ..., quote = FALSE){
    cat("Out-of-sample sensitivity and specificity at outcome thresholds:\n")
    df_ <- data.frame(threshold = x$CVroc$thresholds,
                      sensitivity = x$CVroc$sensitivities,
                      specificity = x$CVroc$specificities)
    print(df_)
}


## Function print.simple_screenr
##
#' An S3 Print Method for \code{screenr} Objects
#'
#' @description \code{print.simple_screenr} is print method for
#' \code{simple_screenr} objects.
#'
#' @param x an object of class \code{simple_screenr}.
#'
#' @param ... optional arguments to \code{print} methods.
#'
#' @return Nothing. Thresholds, specificities and sensitivities are printed as a
#' side effect.
#'
#' @examples
#' data(unicorns)
#' toosimple <- simple_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6,
#'                           data = unicorns)
#' print(toosimple)
#'
#' @export
print.simple_screenr <- function(x, ...){
    cat("\nIn-sample (overly optimistic) sensitivity and specificity:\n")
    df_ <- pROC::coords(x$ISroc, transpose = FALSE)
    df_["threshold"] <- df_["threshold"] + 0.5
    print(df_)
}
