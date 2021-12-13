#################################################################################
##     R Script: helperFunctions.R
##
##      Package: screenr
##
##  Description: Miscellaneous helper functions
##
##       Author: Steve Gutreuter
##               sgutreuter@gmail.com
#################################################################################



## Function keepfirst
##
#' Return Data Frame Rows Having Unique Values in Selected Columns
#'
#' @description \code{keepfirst} extracts those rows of a data frame which have
#' unique values in selected columns.
#'
#' @param x character-valued column name along which the dataframe is sorted.
#'
#' @param colnames a character vector of column names  to identify uniqueness.
#'
#' @param data a data frame.
#'
#' @details
#' The dataframe \code{data} is sorted, and then only those rows which are unique
#' with respect to the values of selected columns.
#'
#' @return A data frame consisting of the rows of \code{data} which are
#' unique with respect to \code{colnames}
keepfirst <- function(x, colnames, data = NULL){
    if(!("data.frame" %in% class(data))) stop("data argument must be a dataframe" )
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


## Function inverseLink
##
#' Compute the Inverses of Binomial Link Functions
#'
#' @description \code{inverseLink} returns the inverse of logit, cloglog and
#' probit link functions for a linear predictor
#'
#' @param lp numeric vector containing the estimated link.
#'
#' @param link (character) name of the link function (one of \verb{"logit"},
#' \verb{"cloglog"} or \verb{"probit"}).
#'
#' @return A numeric vector containing the inverse of the link function for the
#' linear predictor.
#'
#' @details
#' \code{inverseLink} returns the inverses of logit, cloglog and probit link
#' functions, and is provided as a (laborious) way to compute predicted values
#' from the \verb{ModelFit} component of \code{logreg_screenr}-class objects.
#' The \code{predict} methods are a better way to obtain predicted values.
#'
#' @seealso \code{\link[screenr]{predict.logreg_screenr}}
#'
#' @note \code{inverseLink} may not be included in future versions of the \code{screenr}
#' package.
#'
#' @examples
#' ## Make predictions of probability of infection from new observations
#' attach(uniobj2)
#' new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
#'                         testresult = c(NA, NA), Q1 = c(0, 0), Q2 = c(0, 0),
#'                         Q3 = c(0, 0), Q4 = c(0, 0), Q5 = c(0, 1), Q6 = c(0, 1 ),
#'                         Q7 = c(0, 1))
#' mfit <- get_what(from = uniobj2 , what = "ModelFit")
#' coefs <- mfit$coefficients
#' lp <- as.matrix(cbind(rep(1, nrow(new_corns)), new_corns[, 3:9])) %*%
#'            as.matrix(coefs, ncol =  1)
#' (preds <- inverseLink(lp, link = "logit"))
#' ## Note that only the predicted values are returned.
#' @importFrom stats pnorm
#' @export
inverseLink <- function(lp = NULL, link =  NULL){
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



## Function nnt_
##
#' Compute the Ratio of Total Tests Performed  Per Postive Result
#'
#' @description \code{nnt_} computes the anticipated average number of tests
#' performed in order to observe a positive test result.
#'
#' @param dframe a data frame containing columns \code{sensitivity},
#' \code{specificity} and \code{prev}.
#'
#' @importFrom dplyr between
nnt_ <- function(dframe) {
    se <- dframe$sensitivity
    sp <- dframe$specificity
    pv <- dframe$prev
    if(!all(c(dplyr::between(sp, 0, 1 ), dplyr::between(se, 0, 1 ),
              dplyr::between(pv, 0, 1 )))) {
        stop("sensitivity, specificity and prev must be between 0 and 1")
    }
    Etpp <- ((se * pv) + (1 - sp) * (1 - pv)) / (se * pv)
    Epu <- ((1 - se) * pv) / ((pv * (1 - se)) + (sp * (1 - pv) ))
    result <- data.frame(cbind(dframe, Etpp, Epu))
    names(result)[ncol(dframe) + c(1, 2)] <- c("ntpp", "prev_untested")
    result
}



## Function rescale_to_int
##
#' Rescale Positive Vectors or Matrices to Integers
#'
#' @description \code{rescale_to_int} rescales the non-zero elements of
#' real-valued numeric vectors or matrices to integers in the closed
#' interval [1, \verb{max}]. Any zero-valued elements are left unchanged.
#'
#' @param x numeric matrix or vector of non-negative real numbers.
#'
#' @param max the value of largest element in the rescaled integer-valued
#' vector.
#'
#' @param colwise (logical) rescale the matrix by column if \verb{TRUE} (the default) or
#' by row if \verb{FALSE}.
#'
#' @return A matrix of integers corresponding to \code{x} in which smallest
#' \emph{non-zero} element in each column/row is 1 and the largest element is \code{max}. Any
#' elements having value zero are unchanged. If \code{x} is a vector then the result is
#' a \emph{r} x 1 matrix, where \emph{r} is the number of elements in \code{x}.  Otherwise
#' the result is a \emph{r} x \emph{c} matrix where \emph{c} is the number of columns in \code{x}.
#'
#' @seealso \code{\link[scales]{rescale}}
#'
#' @examples
#' x <- c(0.55, 1.21, 0.94, 0, 0.13)
#' rescale_to_int(x, max = 5)
#' @importFrom scales rescale
#' @export
rescale_to_int <- function(x, max, colwise = TRUE){
    if(any(x < 0) | max <= 0) stop("Elements of x must be non-negative")
    if(!(is.vector(x) | is.matrix(x))) stop("x must be a vector or matrix")
    if(is.matrix(x)) {
        x_ <- x
    } else {
        x_ <- as.matrix(x, ncol = 1 )
    }
    i0 <- which(x_ == 0)
    min <- x_[which.min(replace(x_, x_ == 0, NA))]
    x_[i0] <- min
    d_ <- ifelse(colwise, 2, 1)
    y <- apply(x_, d_, function(x_) round(scales::rescale(x_, to = c(1, max))))
    y[i0] <- 0
    y
}


## Function roc_ci
##
#' Compute Bootstrap Confidence Limits for Sensitivities and Specificities
#'
#' @description \code{roc_ci} computes bootstrap confidence intervals from
#' objects of class \code{roc}, as produced by the \code{pROC} package.
#' \code{roc_ci} is simply a convenience wrapper for
#' \code{pROC::ci.thresholds} re-formatted for \code{screenr}.
#'
#' @param object an object of class \code{roc}.
#'
#' @param bootreps number of bootstrap replicates (default = 2000).
#'
#' @param conf.level confidence level for uncertainty intervals
#' (default = 0.95).
#'
#' @param progress type of progress display
#' (see \code{help(pROC::ci.thresholds)}).
#'
#' @param thresholds type of thresholds (see \code{help(pROC::ci.thresholds)}).
#'
#' @param se.min minimum value of sensitivity returned.
#'
#' @return a data frame containing thresholds with sensititives, specificities
#' and uncertainy intervals.
#' @export
roc_ci <- function(object, bootreps = 2000, conf.level = 0.95,
                   progress = "none", thresholds = "local maximas",
                   se.min = 0.5) {
    if(!class(object) == "roc") stop("class(object) must be 'roc'" )
    ci_ <- pROC::ci.thresholds(object,
                               boot.n = bootreps,
                               progress = progress,
                               conf.level = conf.level,
                               thresholds = thresholds)
    res <- cbind(ci_$sensitivity, ci_$specificity)
    threshold <- as.numeric(rownames(res))
    res <- data.frame(cbind(threshold, res))
    rownames(res) <- 1:dim(res)[1]
    names(res) <- c("Threshold", "se.lcl", "Sensitivity",
                    "se.ucl", "sp.lcl", "Specificity", "sp.ucl")
    res <- res[res$Sensitivity >= se.min, c(1, 3, 2, 4, 6, 5, 7)]
    res[is.infinite(res[,1]), 1] <- 0
    res
}


## Function sens_spec
##
#' Compute Sensitivity and Specificity from a 2 x 2 Table
#'
#' @description \code{sens_spec} computes sensitivity and specificity from a
#' 2 x 2 table.
#'
#' @param x a 2 x 2 table, with columns representing frequencies of
#' gold-standard status and rows representing frequencies of status ascertained
#' from testing.  The first row contains frequencies of negative test results
#' and the first column contain frequencies of true negatives.
#'
#' @return a list containing components sensitivity and specificity.
#' Sensitivities and specificities are displayed as proportions rather than
#' percentages.
#'
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


################################   END of FILE   ################################