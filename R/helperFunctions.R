#################################################################################
##       R PROGRAM: helperFunctions.R
##
##         PROJECT: R functions for HIV screening tool development
##
##      WRITTEN BY: Steve Gutreuter
##                  E-mail:  sgutreuter@gmail.com
##
#################################################################################


## Generic function getWhat
##
#' \code{getWhat} is an S3 generic function to extract components of objects
#' produced by functions in the \code{screenr} package
#'
#' @usage \code{getWhat(object, ...)}
#'
#' @seealso \code{link[screenr]{getWhat.glmpathScreenr}}
#' @seealso \code{link[screenr]{getWhat.logisticScreenr}}
#' @export
getWhat <- function(from, ...) {
    UseMethod("getWhat", from)
}


## Generic function ntpp
##
#' \code{ntpp} is an S3 generic function that computes the anticipated number of
#' tests per positive test result and prevalence among the subjects who would be
#' screened out of testing)
#'
#' @usage \code{getWhat(object, ...)}
#
#' @seealso \code{link[screenr]{ntpp.glmpathScreenr}}
#' @seealso \code{link[screenr]{ntpp.logisticScreenr}}
#' @seealso \code{link[screenr]{ntpp.data.frame}}
#'
#' @export
ntpp <- function(object, ... ) {
    UseMethod("ntpp", object)
}


## Function sens_spec
##
#' \code{sens_spec} Computes sensitivity and specificity of a test from a
#' 2 x 2 table
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


## Function inverseLink
##
#' \code{inverseLink} computes inverse of logistic regression link functions
#'
#' Returns the inverse of logit, cloglog and probit link functions for a linear
#' predictor
#'
#' @param lp numeric vector containing the estimated link.
#'
#' @param link (character) name of the link function (one of \verb{"logit"}, \verb{"cloglog"}
#' or \verb{"probit"}).
#'
#' @return A numeric vector containing the inverse of the link function for the
#' linear predictor.
#'
#' @details
#' \code{inverseLink} returns the inverses of logit, cloglog and probit link functions,
#' and is provided as a (laborious) way to compute predicted values from the \verb{ModelFit}
#' component of \code{logisticScreenr}-class objects.  The \code{predict} methods are
#' a better way to obtain predicted values.
#'
#' @seealso \code{\link[screenr]{predict.logisticScreenr}}
#'
#' @note \code{inverseLink} may not be included in future versions of the \code{screenr}
#' package.
#'
#' @examples
#' ## Make predictions of probability of infection from new observations
#' load(bsobj2)
#' class(bsobj2)
#' new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
#'                         testresult = c(NA, NA), Q1 = c(0, 0), Q2 = c(0, 0),
#'                         Q3 = c(0, 0), Q4 = c(0, 0), Q5 = c(0, 1), Q6 = c(0, 1 ))
#' mfit <- getWhat(what = "ModelFit", from = bsobj2 )
#' coefs <- mfit$coefficients
#' lp <- as.matrix(cbind(rep(1, nrow(new_corns)), new_corns[, 3:8])) %*%
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


## Function ntpp.data.frame
##
#' \code{ntpp.data.frame} is a method for computation of the anticipated
#' number of tests per positive test result
#'
#' @param dframe a one-row dataframe containing columns \code{sensitivity},
#' \code{specificity} and \code{prev}.
#'
#' @return a data frame containing the following columns:
#' \describe{
#' \item{\code{sensitivity}}{the sensitivity (proportion)}
#' \item{\code{specificity}}{the specificity (proportion)}
#' \item{\code{prev}}{prevalence proportion of the test condition}
#' \item{\code{ntpp}}{anticipated total tests required per positive result}
#' \item{\code{prev_untested}}{anticipated prevalence proportion among the untested}
#'}
#' @export
ntpp.data.frame <- function(dframe){
    if(!is.data.frame(dframe)) stop("dframe not a data frame")
    if(!"sensitivity" %in% names(dframe))
        stop("dframe does not include sensitivity")
    if(!"specificity" %in% names(dframe))
        stop("dframe does not include specificity")
    if(!"prev" %in% names(dframe))
        stop("dframe does not include prev")
    result <- nnt_(dframe)
    result
}


## Function nnt_
##
#' \code{nnt_} computes the anticpated number of tests per postive from a
#' structured dataframe
#'
#' @param dframe a data frame containing columns \code{sensitivity},
#' \code{specificity} and \code{prev}.
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


## Function keepfirst
##
#' \code{keepfirst} returns only those rows of a dataframe having unique values
#' in selected columns.
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


## Function rescale_to_int
##
#' \code{rescale_to_int} rescales all non-zero elements of a non-negative
#' numeric vector to integers from 1 to \verb{max}
#'
#' @param x numeric vector of positive real numbers.
#'
#' @param max the value of largest element in the rescaled integer-valued
#' vector.
#'
#' @return A vector of integers corresponding to \code{x} in which smallest
#' \emph{non-zero} element is 1 and the largest element is \code{max}. Any
#' elements having value zero are unchanged.
#'
#' @note Any values of 0 in \code{x} are not rescaled, and are preserved in
#' the result.
#'
#' @seealso \code{\link[scales]{rescale}}
#'
#' @examples
#' x <- c(0.5, 1.2, 0.9, 0, 0.1)
#' rescale_to_int(x, max = 5)
#' @import scales
#' @export
rescale_to_int <- function(x, max){
    if(any(x < 0) | max <= 0) stop("Elements of x must be non-negative" )
    x_ <- x
    i0 <- which(x_ == 0)
    min <- x_[which.min(replace(x_, x_ == 0, NA))]
    x_[i0] <- min
    y <- round(scales::rescale(x_, to = c(1, max)))
    y[i0] <- 0
    y
}

################################   END of FILE   ################################
