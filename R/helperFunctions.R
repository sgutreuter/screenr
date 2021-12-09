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
#' attach(uniobj2)
#' new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
#'                         testresult = c(NA, NA), Q1 = c(0, 0), Q2 = c(0, 0),
#'                         Q3 = c(0, 0), Q4 = c(0, 0), Q5 = c(0, 1), Q6 = c(0, 1 ),
#'                         Q7 = c(0, 1))
#' mfit <- getWhat(from = uniobj2 , what = "ModelFit")
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
#' numeric matrix or vector to integers from 1 to \verb{max}
#'
#' @param x numeric matrix or vector of positive real numbers.
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
#' the result is a \emph{r} x \emph{c} matrix where \emph{c} is the number of columns in \code{x}
#'
#' @note Any values of 0 in \code{x} are not rescaled, and are preserved in
#' the result.
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

################################   END of FILE   ################################
