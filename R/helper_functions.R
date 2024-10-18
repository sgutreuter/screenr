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


## Function inverse_link
##
#' Compute the Inverses of Binomial Link Functions
#'
#' @description \code{inverse_link} returns the inverse of logit, cloglog and
#' probit link functions for a linear predictor
#'
#' @param lp numeric vector containing the estimated link.
#'
#' @param link (character) name of the link function (one of \verb{"logit"},
#' \verb{"cloglog"} or \verb{"probit"}).
#'
#' @return \code{inverse_link} returns a numeric vector containing the inverse
#' of the link function for the linear predictor.
#'
#' @details
#' \code{inverse_link} returns the inverses of logit, cloglog and probit link
#' functions, and is provided as a (laborious) way to compute predicted values
#' from the \verb{ModelFit} component of \code{logreg_screenr}-class objects.
#' The \code{predict} methods are a better way to obtain predicted values.
#'
#' @seealso \code{\link[screenr]{predict.logreg_screenr}}
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
#' lp <- as.matrix(cbind(rep(1, nrow(new_corns)), new_corns[, 3:8])) %*%
#'            as.matrix(coefs, ncol =  1)
#' (preds <- inverse_link(lp, link = "logit"))
#' ## Note that only the predicted values are returned.
#' @importFrom stats pnorm
#' @export
inverse_link <- function(lp = NULL, link = c("logit", "cloglog", "probit")){
    link <- match.arg(link)
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
#' @description \code{nnt_} computes the anticipated average number of
#' tests performed in order to observe a positive test result.
#'
#' @param dframe a dataframe containing columns \code{sensitivities},
#' \code{specificities} and \code{prev}.
#'
#' @return \code{nnt_} returns adataframe containing sensitivity,
#' specificity, the anticipated
#' average number of tests required to observe a single positive test
#' result \code{ntpp}, and the prevalence among those screened out
#' of testing \code{pre_untested}.
#'
#' @importFrom dplyr between
nnt_ <- function(dframe) {
    se <- dframe$sensitivities
    sp <- dframe$specificities
    pv <- dframe$prev
    if(!all(c(dplyr::between(sp, 0, 1 ), dplyr::between(se, 0, 1 ),
              dplyr::between(pv, 0, 1 )))) {
        stop("sensitivities, specificities and prev must be between 0 and 1")
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
#' @description \code{rescale_to_int} rescales the \emph{non-zero} elements of
#' real-valued numeric vectors or matrices to integers in the closed
#' interval [1, \verb{max}]. Any zero-valued elements are left unchanged.
#'
#' @param x numeric matrix or vector of non-negative real numbers.
#'
#' @param max the value of largest element in the rescaled integer-valued
#' vector.
#'
#' @param colwise (logical) rescale the matrix by column if \verb{TRUE}
#' (the default) or by row if \verb{FALSE}.
#'
#' @return \code{rescale_to_int} returns a matrix of integers corresponding
#' to \code{x} in which smallest \emph{non-zero} element in each column/row
#' is 1 and the largest element is \code{max}. Any elements having value zero
#' are unchanged. If \code{x} is a vector then the result is an \emph{r} x 1
#' matrix, where \emph{r} is the number of elements in \code{x}.  Otherwise
#' the result is a \emph{r} x \emph{c} matrix where \emph{c} is the number of
#' columns in \code{x}.
#'
#' @seealso \code{\link[scales]{rescale}}
#'
#' @examples
#' x <- c(0.55, 1.21, 0.94, 0, 0.13)
#' rescale_to_int(x, max = 5)
#' @importFrom scales rescale
#' @export
rescale_to_int <- function(x, max, colwise = TRUE){
    if(!(is.vector(x) | is.matrix(x))) stop("x must be a vector or matrix")
    if(any(x < 0) | max <= 0) stop("Elements of x must be non-negative")
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
#' @param bootreps number of bootstrap replicates.  Default: 4000.
#'
#' @param conf_level confidence level for uncertainty intervals.
#' Default: 0.95.
#'
#' @param progress character-valued type of progress display
#' (see \code{help(pROC::ci.thresholds)}). Default \verb{"none"}.
#'
#' @param thresholds type of thresholds (see \code{help(pROC::ci.thresholds)}).
#'
#' @param se_min minimum value of sensitivity returned. Default: 0.8.
#'
#' @param ... additional arguments passed to \code{pROC::ci.thresholds}
#'
#' @return \code{roc_ci} returns a dataframe containing thresholds with
#' their sensititives, specificities and uncertainy intervals.
#'
#' @seealso \code{\link[pROC]{ci.thresholds}}
#'
#' @export
roc_ci <- function(object, ..., bootreps = 4000, conf_level = 0.95,
                   progress = "none", thresholds = "local maximas",
                   se_min = 0.8) {
    if(!inherits(object, "roc")) stop("class(object) must be 'roc'" )
    ci_ <- pROC::ci.thresholds(object,
                               boot.n = bootreps,
                               progress = progress,
                               conf.level = conf_level,
                               thresholds = thresholds)
    res <- data.frame(cbind(ci_$sensitivity, ci_$specificity))
    thresholds <- as.numeric(rownames(res))
    res <- data.frame(thresholds = thresholds, res)
    rownames(res) <- 1:dim(res)[1]
    names(res) <- c("thresholds", "se.lcl", "sensitivities",
                    "se.ucl", "sp.lcl", "specificities", "sp.ucl")
    res <- res[res$sensitivities >= se_min, c(1, 3, 2, 4, 6, 5, 7)]
    res[is.infinite(res[,1]), 1] <- 0
    res
}


## Function se_sp_max
##
#' Return a Simplified Dataframe of Sensitivity and Specificity
#'
#' @description Given a dataframe containing multiple values of specificity for
#' each value of sensitivity, return only the rows containing the largest value of
#' specificity for each unique value of sensitivity.
#'
#' @param object a dataframe containing at least columns named \code{sensitivities}
#' and \code{specificities}
#'
#' @return \code{se_sp_max} returns a dataframe which is a subset of \code{object}
#' containing only those rows for which specificity was the maximum for each
#' unique value of sensitivity.
#'
#' @import magrittr
#' @importFrom dplyr group_by summarise arrange desc right_join distinct
se_sp_max <- function(object) {
    stopifnot(is.data.frame(object))
    stopifnot(all(c("sensitivities", "specificities") %in% names(object)))
    sensitivities <- specificities <- NULL
    object <- dplyr::distinct(object)
    mxse <- object %>%
        dplyr::group_by(sensitivities) %>%
        dplyr::summarise(specificities = max(specificities)) %>%
        dplyr::arrange(dplyr::desc(sensitivities))
    dplyr::right_join(object, mxse)
}


## Function sens_spec_plus
##
#' Compute Sensitivity, Specificity and a Few Friends
#'
#' @description \code{sens_spec_plus} computes sensitivity, specificity and a
#' few friends from a gold standard and testing results. \code{sens_spec_plus}
#' is a convenience wrapper for \code{epiR::epi.tests}.
#'
#' @param test numeric vector containing testing results, coded as 0 for
#' negative and 1 for positive.
#'
#' @param gold numeric vector containing the gold standard, coded as 0 for
#' negative and 1 for positive.
#'
#' @param method type of uncertainty interval
#' (\verb{"exact", "wilson", "agresti", "clopper-pearson" or"jeffreys"}). Default:
#' \verb{"exact"}.
#'
#' @param conf_level confidence level, a numeric value between 0 and 1.
#' Default: 0.95.
#'
#' @return \code{sens_spec_plus} returns a list containing components
#' \code{table} and \code{ests}:
#' \describe{
#' \item{\code{table}}{a 2 x 2 table which is the anti-transpose of the
#' result produced by \code{base::table(gold, test)}.}
#' \item{\code{ests}}{a dataframe containing the apparent (test-based) and true
#' positive proportions, sensitivity, specificity, positive predictive value
#' (PPV), negative predictive value (NPV), and the lower and upper uncertainty
#' limits for each.}
#' }
#'
#' @seealso \code{\link[epiR]{epi.tests}}
#'
#' @examples
#' Gold <- rbinom(20, 1, 0.50)
#' Test <- Gold; Test[c(3, 5, 9, 12, 16)] <- 1 - Test[c(3, 5, 9, 12, 16)]
#' sens_spec_plus(test = Test, gold = Gold, method = "jeffreys")
#'
#' @importFrom epiR epi.tests
#' @export
sens_spec_plus <- function(test = NULL, gold = NULL,
                           method = c("exact", "jeffreys", "wilson", "agresti",
                                      "clopper-pearson"),
                           conf_level =  0.95){
    if(!all((test %in% c(0, 1)) | is.na(test)))
        stop("Values of test must be 0, 1 or NA")
    if(!all((gold %in% c(0, 1)) | is.na(gold)))
        stop("Values of gold must be 0, 1 or NA")
    if(!(conf_level > 0 & conf_level < 1 )) stop("conf.level not in (0,1)")
    cimethod <- match.arg(method)
    .tst <- ordered(test, levels = c(1, 0), labels = c("Pos", "Neg"))
    .gld <- ordered(gold, levels = c(1, 0), labels = c("Pos", "Neg"))
    .tbl <- table(.tst, .gld)
    .res <- epiR::epi.tests(.tbl, method = cimethod, conf.level = conf_level)
    .table <- .res$tab
    colnames(.table) <- c("      True +", "      True -", "     Total")
    .ests <- .res$detail[c(1:4, 9:10), c(2:4)]
    rownames(.ests) <- c("Apparent_positivity", "True_positivity", "Sensitivity",
                        "Specificity", "PPV", "NPV")
    result <- list(table = .table, ests = .ests)
    attr(result, "Interval_type") <- method
    attr(result, "Confidence_level") <- conf_level
    result
}
################################   END of FILE   ################################
