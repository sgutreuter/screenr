################################################################################
##     R Script: ntpp.R
##
##      Package: screenr
##
##  Description: S3 methods for number of tests per positive result
################################################################################


## Generic function ntpp
##
#' \code{ntpp} is an S3 generic function that computes the anticipated number of
#' tests per positive test result and prevalence among the subjects who would be
#' screened out of testing)
#'
#' @param object an object from which to compute the number of tests
#' per test positive test results.
#'
#' @param ... additional arguments.
#'
#' @seealso \code{link[screenr]{ntpp.lasso_screenr}}
#' @seealso \code{link[screenr]{ntpp.logreg_screenr}}
#' @seealso \code{link[screenr]{ntpp.data.frame}}
#' @seealso \code{link[screenr]{ntpp.simple_screenr}}
#'
#' @export
ntpp <- function(object, ...) {
    UseMethod("ntpp")
}


## Function ntpp.easy_tool
##
#' \code{ntpp.easy_tool} is a method for computation of the anticipated
#' number of tests per positive test result
#'
#' @param object an \code{easy_tool}-class object produced by \code{easy_tool}.
#'
#' @param prev an optional prevalence proportion for the test outcome; if missing
#' the prevalence is obtained from \code{object}.
#'
#' @param ... optional arguments to \code{ntpp} methods.
#'
#' @return A data frame containing the following columns:
#' \describe{
#' \item{\verb{sensitivity}}{The sensitivity (proportion) of the screener.}
#' \item{\verb{specificity}}{The specificity (proportion) of the screener.}
#' \item{\verb{ntpp}}{the number of tests required to discover
#' a single positive test result.}
#' \item{\verb{prev_untested}}{The prevalence proportion of the test
#' condition among those who are screened out of testing.}
#' }
#'
#' @details
#' The anticipated number of tests needed to observe a single positive test
#' result is a function of sensitivity, specificity and the prevalence proportion
#' of the condition being tested.
#'
#' @examples
#' attach(uniobj1)
#' tool <- easy_tool(uniobj1, max = 3, crossval = TRUE)
#' ntpp(tool)
#'
#' @export
ntpp.easy_tool <- function(object, ..., prev = NULL) {
     if(!class(object) == "easy_tool")
         stop("object not of class easy_tool")
     if(is.null(prev)) prev <- mean(object$Scores$response, na.rm = TRUE)
     ssp <- data.frame(sensitivity = object$ROC$sensitivities,
                       specificity = object$ROC$specificities)
     ssp <- cbind(ssp, rep(prev, dim(ssp)[1]))
     names(ssp) <- c("sensitivity", "specificity", "prev")
     result <- nnt_(ssp)
     result
}


## Function ntpp.lasso_screenr
##
#' \code{ntpp.lasso_screenr} is a method for computation of the anticipated
#' number of tests per positive test result
#'
#' @param object a \code{lasso_screenr}-class object produced by \code{lasso_screenr}.
#'
#' @param model (character) select the model which produced the
#' minimum AIC (\verb{"minAIC"}, the default) or minimum BIC (\verb{"minBIC"}).
#'
#' @param type (character) one of \verb{"cvResults"} (the default) or
#' \verb{"isResults"} to specify \emph{k}-fold cross-validated or in-sample
#' receiver-operating characteristics, respectively.
#'
#' @param prev an optional prevalence proportion for the test outcome; if missing
#' the prevalence is obtained from \code{object}.
#'
#' @param ... optional arguments to \code{ntpp} methods.
#'
#' @return A data frame containing the following columns:
#' \describe{
#' \item{\verb{sensitivity}}{The sensitivity (proportion) of the screener.}
#' \item{\verb{specificity}}{The specificity (proportion) of the screener.}
#' \item{\verb{ntpp}}{the number of tests required to discover
#' a single positive test result.}
#' \item{\verb{prev_untested}}{The prevalence proportion of the test
#' condition among those who are screened out of testing.}
#' }
#'
#' @details
#' The anticipated number of tests needed to observe a single positive test
#' result is a function of sensitivity, specificity and the prevalence proportion
#' of the condition being tested.
#'
#' @examples
#' attach(uniobj1)
#' ntpp(uniobj1)
#'
#' @export
ntpp.lasso_screenr <- function(object, ..., model = "minAIC", type = "cvResults",
                                prev = NULL) {
     if(!class(object) == "lasso_screenr")
         stop("object not of class lasso_screenr")
     if(!type %in% c("cvResults", "isResults"))
         stop("type must be 'cvResults' or 'isResults'" )
     if(is.null(prev)) prev <- object$Prevalence
     ssp <- data.frame(sensitivity = object[[type]][[model]][["ROC"]][["sensitivities"]],
                       specificity = object[[type]][[model]][["ROC"]][["specificities"]])
     ssp <- cbind(ssp, rep(prev, dim(ssp)[1]))
     names(ssp) <- c("sensitivity", "specificity", "prev")
     result <- nnt_(ssp)
     result
}


## Function ntpp.logreg_screenr
##
#' \code{ntpp.logreg_screenr} is a method for computation of the anticipated
#' number of tests per positive test result
#'
#' @param object a \code{logreg_screenr}-class object produced by \code{logreg_screenr}.
#'
#' @param type (character) one of \verb{"cvResults"} (the default) or
#' \verb{"isResults"} to specify \emph{k}-fold cross-validated or in-sample
#' receiver-operating characteristics, respectively.
#'
#' @param prev an optional prevalence proportion for the test outcome; if missing
#' the prevalence is obtained from \code{object}.
#'
#' @param ... optional arguments to \code{ntpp} methods.
#'
#' @return A data frame containing the following columns:
#' \describe{
#' \item{\verb{sensitivity}}{The sensitivity (proportion) of the screener.}
#' \item{\verb{specificity}}{The specificity (proportion) of the screener.}
#' \item{\verb{ntpp}}{the number of tests required to discover
#' a single positive test result.}
#' \item{\verb{prev_untested}}{The prevalence proportion of the test
#' condition among those who are screened out of testing.}
#' }
#'
#' @details
#' The anticipated number of tests needed to observe a single positive test
#' result is a function of sensitivity, specificity and the prevalence proportion
#' of the condition being tested.
#'
#' @examples
#' attach(uniobj2)
#' ntpp(uniobj2)
#'
#' @export
ntpp.logreg_screenr <- function(object, ..., type = "cvResults",
                                prev = NULL) {
     if(!class(object) == "logreg_screenr")
         stop("object not of class logreg_screenr")
     if(!type %in% c("cvResults", "isResults"))
         stop("type must be 'cvResults' or 'isResults'" )
     if(is.null(prev )) prev <- object$Prevalence
     if(type == "cvResults") {
         x <- "CVroc"
     } else {
         x <- "ISroc"
     }
     ssp <- data.frame(sensitivity = object[[x]][["sensitivities"]],
                       specificity = object[[x]][["specificities"]],
                       prev = prev)
     result <- nnt_(ssp)
     result
}


## Function ntpp.data.frame
##
#' \code{ntpp.data.frame} is a method for computation of the anticipated
#' number of tests per positive test result
#'
#' @param object a one-row dataframe containing columns \code{sensitivity},
#' \code{specificity} and \code{prev}.
#'
#' @param ... optional arguments to \code{ntpp} methods.
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
ntpp.data.frame <- function(object, ...){
    if(!is.data.frame(object)) stop("object not a data frame")
    if(!"sensitivity" %in% names(object))
        stop("object does not include sensitivity")
    if(!"specificity" %in% names(object))
        stop("object does not include specificity")
    if(!"prev" %in% names(object))
        stop("object does not include prev")
    result <- nnt_(object)
    result
}


## Function ntpp.simple_screenr
##
#' \code{ntpp.simple_screenr} is a method for computation of the anticipated
#' number of tests per positive test result
#'
#' @param object a \code{simple_screenr}-class object produced by \code{simple_screenr}.
#'
#' @param prev an optional prevalence proportion for the test outcome; if missing
#' the prevalence is obtained from \code{object}.
#'
#' @param ... optional arguments to \code{ntpp} methods.
#'
#' @return A data frame containing the following columns:
#' \describe{
#' \item{\verb{sensitivity}}{The sensitivity (proportion) of the screener.}
#' \item{\verb{specificity}}{The specificity (proportion) of the screener.}
#' \item{\verb{ntpp}}{the number of tests required to discover
#' a single positive test result.}
#' \item{\verb{prev_untested}}{The prevalence proportion of the test
#' condition among those who are screened out of testing.}
#' }
#'
#' @details
#' The anticipated number of tests needed to observe a single positive test
#' result is a function of sensitivity, specificity and the prevalence proportion
#' of the condition being tested.
#'
#' @examples
#' data(unicorns)
#' toosimple <- simple_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
#'                            data = unicorns)
#' ntpp(toosimple)
#'
#' @export
ntpp.simple_screenr <- function(object, ..., prev = NULL) {
     if(!class(object) == "simple_screenr")
         stop("object not of class simple_screenr")
     if(is.null(prev)) prev <- object$Prevalence
     ssp <- data.frame(sensitivity = object[["ISroc"]][["sensitivities"]],
                       specificity = object[["ISroc"]][["specificities"]])
     ssp <- cbind(ssp, rep(prev, dim(ssp)[1]))
     names(ssp) <- c("sensitivity", "specificity", "prev")
     result <- nnt_(ssp)
     result
}
