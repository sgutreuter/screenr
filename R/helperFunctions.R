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
#' Sensitivities and specificities are displayed as proportions rather than
#' percentages.
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
#' @param link Character link function (one of \verb{"logit"}, \verb{"cloglog"}
#' or \verb{"probit"})
#' @param lp Numeric vector containin the linear predictor
#'
#' @return A numeric vector containing the inverse of the link function for the
#' linear predictor.
#' @importFrom stats pnorm
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
#' @return A data frame containing threshold scores, sensitivities and
#' specificities. Sensitivities and specificities are displayed as proportions
#' rather than percentages.
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
    }
    res <- data.frame(threshold = th,
                      sensitivity = obj$sensitivities,
                      specificity = obj$specificities)
    res
}


#' Expected Number of Tests Required per Positive Test Result
#'
#' Compute the expected number of tests which need to be performed in order
#' to identify the first positive test result and the expected prevalence
#' among the untested given implementation of test screening options having
#' the specified values of sensitivity and specificity.
#'
#' @param x A data frame containing columns "sensitivity" and
#' "specificity", or an object of class 'simplescreenr' or 'binomscreenr'.
#' @param prev Numeric proportion of the population expressing positive test
#' results.  \code{prev} is \emph{optional} for class 'simplescreenr' and
#' 'binomscreenr' objects, for which the default is the prevalence of the test
#' condition in the training sample.
#'
#' @return A data frame containing the following columns:
#' \describe{
#' \item{\verb{sensitivity}}{The sensitivity (proportion) of the screener.}
#' \item{\verb{specificity}}{The specificity (proportion) of the screener.}
#' \item{\verb{E(Tests/Pos)}}{The expected number of tests required to discover
#' a single positive test result.}
#' \item{\verb{E(prev|untested)}}{The expected prevalence proportion of the test
#' condition among those who are screened out of testing.}
#' }
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
    Epu <- ((1 - ss[["sensitivity"]]) * prev) / ((1 - ss[["sensitivity"]]) * prev
        + ss[["specificity"]] * (1 - prev))
    result <- data.frame(cbind(ss, Etpp, Epu))
    names(result)[ncol(ss) + c(1, 2)] <- c("E(Tests/Pos)", "E(prev|untested)")
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
