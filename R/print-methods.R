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
#' A Print Method for \code{easy_tool}-Class Objects
#'
#' @description \code{print.easy_tool} is a print method.
#'
#' @param x an object of class \code{easy_tool}.
#'
#' @param ... optional arguments to \code{print} methods.
#'
#' @seealso \code{get_what.easy_tool(from, what)} for \code{what =} \verb{"ROCci"}.
#'
#' @examples
#' attach(uniobj1)
#' print(uniobj1)
#' @export
print.easy_tool <- function(x, ...){
    if(!("easy_tool" %in% class(x)))
        stop("x not easy_tool class")
    cat("Function call:\n")
    print(x$Call)
    cat("\nType:\n")
    print(x$Type )
    cat("\nQuestion weights:\n")
    print(x$QuestionWeights )
}


## Function print.lasso_screenr
##
#' A Print Method for \code{lasso_screenr}-Class Objects
#'
#' @description \code{print.lasso_screenr} is a print method for
#' \code{lasso_screenr}-class objects.
#'
#' @param x an object of class \code{lasso_screenr}
#'
#' @param ... optional arguments to \code{print} methods.
#'
#' @seealso \code{get_what.lasso_screenr(from, what)} for \code{what =} \verb{"ROCci"}.
#'
#' @examples
#' attach(uniobj1)
#' print(uniobj1)
#' @export
print.lasso_screenr <- function(x, ...){
    if(!("lasso_screenr" %in% class(x)))
        stop("x not lasso_screenr class")
    cat("Function call:\n")
    print(x$Call)
    cat("\nglmpath x:\n")
    glmpath::print.glmpath(x$glmpathObj)
    cat("\nTest positivity:", x$Prevalence, "\n")
}


## Function print.logreg_screenr
##
#' A Print Method for \code{logreg_screenr}-Class Objects
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
#' @seealso \code{get_what.logreg_screenr(from, what)} for \code{what =} \verb{"ROCci"}.
#'
#' @examples
#' attach(uniobj2)
#' print(uniobj2)
#'
#' @export
print.logreg_screenr <- function(x, ..., quote = FALSE){
    if(!("logreg_screenr" %in% class(x))) stop("x not logreg_screenr class")
    cat("Out-of-sample sensitivity and specificity at outcome thresholds:\n")
    df_ <- data.frame(threshold = x$CVroc$thresholds,
                      sensitivity = x$CVroc$sensitivities,
                      specificity = x$CVroc$specificities)
    print(df_)
}


## Function print.simple_screenr
##
#' A Print Method for \code{simple_screenr}-Class Objects
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
#' @seealso \code{get_what.simple_screenr(from, what)} for \code{what =} \verb{"ROCci"}.
#'
#' @examples
#' data(unicorns)
#' toosimple <- simple_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6,
#'                           data = unicorns)
#' print(toosimple)
#'
#' @export
print.simple_screenr <- function(x, ...){
    if(!("simple_screenr" %in% class(x))) stop("x not a simple_screenr object")
    cat("\nIn-sample (overly optimistic) sensitivity and specificity:\n")
    df_ <- pROC::coords(x$ISroc, transpose = FALSE)
    df_["threshold"] <- df_["threshold"] + 0.5
    print(df_)
}
