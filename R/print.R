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


## Function print.easyTool
##
#' \code{print.easyTool} is a print method for \code{easyTool} objects
#'
#' @param x an object of class \code{easyTool}
#'
#' @param ... optional arguments to \code{print} methods.
#'
#' @examples
#' attach(uniobj1)
#' print(uniobj1)
#' @export
print.easyTool <- function(x, ...){
    if(!("easyTool" %in% class(x)))
        stop("x not easyTool class")
    cat("Function call:\n")
    print(x$Call)
    cat("\nType:\n")
    print(x$Type )
    cat("\nQuestion weights:\n")
    print(x$QuestionWeights )
}


## Function print.glmpathScreenr
##
#' code{print.glmpathScreenr} is a print method for \code{glmpathScreenr} objects
#'
#' @param x an object of class \code{glmpathScreenr}
#'
#' @param ... optional arguments to \code{print} methods.
#'
#' @examples
#' attach(uniobj1)
#' print(uniobj1)
#' @export
print.glmpathScreenr <- function(x, ...){
    if(!("glmpathScreenr" %in% class(x)))
        stop("x not glmpathScreenr class")
    cat("Function call:\n")
    print(x$Call)
    cat("\nglmpath x:\n")
    glmpath::print.glmpath(x$glmpathObj)
    cat("\nPrevalence:", x$Prevalence, "\n")
}


## Function print.logisticScreenr
##
#' \code{print.logisticScreenr} is an S3 print method for \code{logisticScreenr} objects.
#'
#' @param x an object of class \code{logisticScreenr}.
#'
#' @param quote logical indicator for whether or not strings should be printed.
#'
#' @param ... optional arguments to \code{print} methods.
#'
#' @return
#' Nothing. Thresholds, specificities and sensitivities are printed as a
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
#' @examples
#' attach(uniobj2)
#' print(uniobj2)
#'
#' @export
print.logisticScreenr <- function(x, ..., quote = FALSE){
    if(!("logisticScreenr" %in% class(x))) stop("x not logisticScreenr class")
    cat("Out-of-sample sensitivity and specificity at outcome thresholds:\n")
    df_ <- data.frame(threshold = x$CVroc$thresholds,
                      sensitivity = x$CVroc$sensitivities,
                      specificity = x$CVroc$specificities)
    print(df_)
}


## Function print.simpleScreenr
##
#' \code{print.simpleScreenr} is an S3 print method for \code{simpleScreenr} objects.
#'
#' @param x an object of class \code{simplescreenr}.
#'
#' @param ... optional arguments to \code{print} methods.
#'
#' @return Nothing. Thresholds, specificities and sensitivities are printed as a
#' side effect.
#'
#' @examples
#' data(unicorns)
#' toosimple <- simpleScreenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6,
#'                           data = unicorns)
#' print(toosimple)
#'
#' @export
print.simpleScreenr <- function(x, ...){
    if(!("simpleScreenr" %in% class(x))) stop("x not a simpleScreenr object")
    cat("\nIn-sample (overly optimistic) sensitivity and specificity:\n")
    df_ <- pROC::coords(x$ISroc, transpose = FALSE)
    df_["threshold"] <- df_["threshold"] + 0.5
    print(df_)
}
