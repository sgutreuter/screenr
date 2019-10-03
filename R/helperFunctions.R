#################################################################################
##       R PROGRAM: helperFunctions.R
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


#' Compute Sensitivity and Specificity from a 2 x 2 Table
#'
#' Computes sensitivity and specificity of a test.
#'
#' @param x A 2 x 2 table, with the numbers of negative test results appearing
#' first in both rows and columns.
#' @return A list containing components sensitivity and specificity.
#' @examples
#' TestResults <- ordered(c(0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0))
#' TrueResults <- ordered(c(0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0))
#' sens_spec(table(TestResults, TrueResults))
#' @export
sens_spec <- function(x){
    if(!class(x) == "table") stop('arg class is not "table"')
    if(!(dim(x)[1] == 2L & dim(x)[2] == 2L)) stop("arg not a 2x2 table")
    spec = x[1, 1] / sum(x[ , 1])
    sens = x[2, 2] / sum(x[ , 2])
    result = list(sensitivity = sens, specificity = spec)
    result
}


#' Compute Inverses of Binomial Regression Link Functions
#'
#' Returns the inverse of logit, cloglog and probit link functions for a linear
#' predictor
#'
#'
#' @param link Character link function (one of "logit", "cloglog" or "probit")
#' @param lp Numeric linear predictor
#'
#' @return The inverse of the link function for the linear predictor.
#' @export
inverseLink <- function(link, lp){
    if(!link %in% c("logit", "cloglog", "probit")) stop("Bad link specification")
    if(link == "logit"){
        p <- exp(lp) / (1 + exp(lp))
    }else{
        if(link == "cloglog"){
            p <- 1 - exp(-exp(lp))
        } else {
            p <- pnorm(lp)
        }
    }
    p
}


#' Extract ROCs from "binomscreenr" or "simplescreenr" Objects
#'
#' Extract the receiver operating characteristics from an object of class
#' "simplescreenr" or "binomscreenr".  This is a convenience function to enable
#' easy use and export of the ROC.
#'
#' @param x An object of class "binomscreenr" or "simplescreenr".
#'
#' @return A data frame containing the threshold scores, sensitivity and
#' specificity.
#'
#' @references
#' Fawcett T. An introduction to ROC analysis. Pattern Recognition Letters. 2006.
#' 27(8):861-874.
#' \url{https://doi.org/10.1016/j.patrec.2005.10.010}
#'
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
#' data(unicorns)
#' unitool <- binomialScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5,
#'                              data = unicorns, Nfolds = 20)
#' (uniROC <- getROC(unitool))
#'
#' @export
getROC <- function(x){
    if(!class(x) %in% c("binomscreenr", "simplescreenr"))
        stop("x not a binomscreenr or simplescreenr object.")
    if(class(x) == "binomscreenr"){
        obj <- x$CVroc
        th <- obj$thresholds
    } else {
        obj <- x$ISroc
        th <- ceiling(obj$thresholds)
>>>>>>> fixROC
    }
    res <- data.frame(threshold = th,
                      sensitivity = obj$sensitivities,
                      specificity = obj$specificities)
    res
}


#' Expected Number of Tests Required per Positive Test Result
#'
#' Compute the expected number of tests which need to be performed in order
#' to identify the first positive test result, and the expected number of
#' false positives among that number of tests.
#'
#' @param x A data frame containing columns "sensitivity" and
#' "specificity", or an object of class 'simplescreenr' or 'binomscreenr'.
#' @param prev Numeric proportion of the population expressing positive test
#' results.  \code{prev} is optional for class 'simplescreenr' and
#' 'binomscreenr' objects, and defaults to prevalence in the training sample
#' if not specified.
#'
#' @return A data frame containing sensitivity, specificity, the expected
#' number of tests required to observe a single positive test result and,
#' among those, the expected number of false negatives per positive test
#' result.
#'
#' @examples
#' data(unicorns)
#' unitool <- binomialScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5,
#'                              data = unicorns, Nfolds = 20)
#' testCounts(unitool)
#'
#' @export
testCounts <- function(x = NULL, prev = NULL){
    if("binomscreenr" %in% class(x)){
        ss <- data.frame(sensitivity = x$CVroc$sensitivities,
                         specificity = x$CVroc$specificities)
        if(is.null(prev)) prev <- x$Prevalence
    } else {
        if("simplescreenr" %in% class(x)){
            ss <- data.frame(sensitivity = x$ISroc$sensitivities,
                             specificity = x$ISroc$specificities)
            if(is.null(prev)) prev <- x$Prevalence
        } else {
            ss <- x
            if(!("sensitivity" %in% names(ss))) stop("No column 'sensitivity")
            if(!("specificity" %in% names(ss))) stop("No column 'specificity")
            if(any(ss[["sensitivity"]] < 0 | ss[["sensitivity"]] > 1))
                stop("sensitivity not in (0,1)")
            if(any(ss[["specificity"]] < 0 | ss[["specificity"]] > 1))
                stop("specificity not in (0,1)")
            if(is.null(prev)) stop("Argument prev missing")
            if(!(prev > 0 & prev < 1)) stop("prev must be in (0,1)")
        }
    }
    Etpp <- ((ss[["sensitivity"]] * prev) +
            (1 - ss[["specificity"]]) * (1 - prev)) / (ss[["sensitivity"]] * prev)
    Efn <- (1 - ss[["sensitivity"]]) * Etpp
    result <- data.frame(cbind(ss, Etpp, Efn))
    names(result)[ncol(ss) + c(1, 2)] <- c("E(Tests/Pos)",
                                           "E(FalseNegs/Pos)")
    result
}


#' Return Data Frame Rows Having Unique Values in Selected Columns
#'
#' Sort a dataframe and then retain those rows which are unique with respect
#' to the values of selected columns.
#' @param x character-valued column name along which the dataframe is sorted.
#' @param colnames a character vector of column names  to identify uniqueness.
#' @param data a data frame.
#'
#' @return A data frame consisting of the rows of \code{data} which are
#' unique with respect to \code{colnames}
keepfirst <- function(x, colnames, data = NULL){
    data <- data[order(data[[x]]), ]
    res <- data[1 ,]
    for(i in 2:(nrow(data) - 1)){
        if(all(data[i, colnames] ==
               data[i + 1, colnames])){
            next
        } else {
            res <- rbind(res, data[i + 1, ])
        }
    }
    res
}


################################   END of FILE   ################################
