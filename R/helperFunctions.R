#################################################################################
##  R CODE COLLECTION: helperFunctions.R
##
##            PACKAGE: screenr
##
##        DESCRIPTION: R functions for HIV screening tool development
##
##         WRITTEN BY: Steve Gutreuter
##                     sgutreuter@gmail.com
##
#################################################################################

## Function sens_spec
##
#' Compute Sensitivity and Specificity from a 2 x 2 Table
#'
#' @param x a 2 x 2 table, with columns representing frequencies of
#' gold-standard status and rows representing frequencies of status ascertained
#' from testing.  The first row contains frequencies of negative test results
#' and the first column contain frequencies of true negatives.
#' @return a list containing components sensitivity and specificity.
#' Sensitivities and specificities are displayed as proportions rather than
#' percentages.
#' @examples
#' Gold <- rbinom(20, 1, 0.50)
#' Test <- Gold; Test[c(3, 9, 12, 16)] <- 1 - Test[c(3, 9, 12, 16)]
#' sens_spec(table(Test, Gold))
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
#' @param lp numeric vector containing the estimated linear predictor.
#' @param link character link function (one of \verb{"logit"}, \verb{"cloglog"}
#' or \verb{"probit"}).
#'
#' @return A numeric vector containing the inverse of the link function for the
#' linear predictor.
#' @importFrom stats pnorm
#' @export
inverseLink <- function(lp, link){
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


#### TODO: Finish inclusion of logisticScreener and simpleScreenr objects.

## Function getWhat
##
#' \code{getWhat} extracts elements from objects of class \code{glmpathScreener} and
#' \code{simpleScreenr}.
#'
#' @param from an object of class glmpathScreener.
#' @param what (character) the element to be extracted; valid values are
#' \verb{"cvPreds"} (cross-validated predicted probabilities),
#' \verb{"cvROC"} (cross validated \code{roc}-class object),
#' \verb{"isPreds"} (predicted probabilities from the training set),
#' \verb{"isROC"} (\code{roc}-class object from the training set) and
#' \verb{"glmpathObj"} (the entire \code{glmpath}-class object).
#' @param model (character) the model from which \code{what} is desired
#' (\verb{"minAIC"} or \verb{"minBIC"}).  The value of \code{model} does not
#' matter for \code{what =} \verb{glmpathObj}, but one of the two valid values
#' must be specified (yes, that is a bit weird).
#' @return the objects specified by \code{what}.
#' @details This function is used to access the specified objects for those who
#' need to perform computations not provided by the methods for class
#' \code{glmpathScreener}.
#' @export
getWhat <- function(from,
                    what = c("cvPreds", "isPreds",
                             "cvROC", "isROC", "glmpathObj"),
                    model = "minAIC"){
    if(!("glmpathScreener" %in% class(from)))
        stop("Object not glmpathScreener class")
    if(!what %in% c("cvPreds", "isPreds", "cvROC", "isROC", "glmpathObj"))
        stop("Invalid what; valid choices are 'cvPreds', 'cvROC, 'isPreds and 'isROC'.")
    if(!what == "glmpathObj"){
        if(!model %in% c("minAIC", "minBIC"))
        stop("Specify 'minAIC' or 'minBIC' for model")
        pfx <- substring(what, 1, 2)
        typ <- substring(what, 3)
        res <- paste0(pfx, "Results")
        x <- from[[res]][[model]][[typ]]
    } else {
        x <- from$glmpathObj
    }
    invisible(x)
}


## Function testCounts
##
#' Expected Number of Tests Required per Positive Test Result
#'
#' Compute the expected number of tests which need to be performed in order
#' to identify the first positive test result and the expected prevalence
#' among the untested given implementation of test screening options having
#' the specified values of sensitivity and specificity.
#'
#' @param x a data frame containing columns \code{sensitivity} and
#' \code{specificity}, or an object of class \code{glmpathScreener} or
#' \code{simplescreenr}.
#' @param prev numeric proportion of the population expressing positive test
#' results.  \code{prev} is \emph{optional} for class \code{glmpathScreener} or
#' \code{simplescreenr} objects, for which the default is the prevalence of
#' the test condition in the training sample.
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
#' unitool <- logisticScreener(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5,
#'                              data = unicorns, Nfolds = 20)
#' testCounts(unitool)
#'
#' @export
testCounts <- function(x = NULL, prev = NULL){
    if("glmpathScreener" %in% class(x)){
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


#' Rescale a strictly positive vector of real numbers to integers
#'
#' A convenience wrapper function which rescales a strictly positive (all
#' elements are greater than zero) vector to a vector of
#' integers ranging from 1 to \code{max}.
#'
#' @param x numeric vector of non-negative real numbers.
#' @param max the value of largest element in the rescaled integer-valued
#' vector.
#'
#' @return A vector of integers corresponding to \code{x} in which smallest
#' element is 1 and the largest element is \code{max}
#'
#' @seealso \code{\link[scales]{rescale}}
#'
#' @import scales
#' @export
rescale_to_Int <- function(x, max) {
    if(any(x <= 0) | max <= 0) stop("x must consist of non-negative numbers" )
    y <- round(scales::rescale(x, to = c(1, max)))
    y
}

################################   END of FILE   ################################
