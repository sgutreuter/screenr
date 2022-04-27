################################################################################
##     R Script: ntpp.R
##
##      Package: screenr
##
##  Description: S3 methods for number of tests per positive result
################################################################################


## Generic function ntpp
##
#' An S3 Method to Compute the Ratio of Total Tests to Positive Results
#'
#' @description \code{ntpp} computes the ratio of the total number of tests
#' performed per positive test result.
#'
#' @param object an object from which to compute the number of tests
#' per test positive test results.
#'
#' @param ... additional arguments.
#'
#' @details
#' The anticipated number of tests required to detect a single positive
#' \emph{nntp} is given by
#' \deqn{nntp = (SeP + (1 - Sp)(1 - P)) / SeP}
#' where \emph{Se} is sensitivity, \emph{P} is prevalence and \emph{Sp} is
#' specificity. The anticipated prevalence among those screened out is given by
#' \deqn{Puntested = ((1 - Se)P) / ((1 - Se)P + Sp (1 - P))}
#'
#' @seealso
#' \code{\link{ntpp.lasso_screenr}}
#' \code{\link{ntpp.logreg_screenr}}
#' \code{\link{ntpp.data.frame}}
#' \code{\link{ntpp.simple_screenr}}
#'
#' @export
ntpp <- function(object, ...) {
    UseMethod("ntpp")
}


## Function ntpp.easy_tool
##
#' Compute the Ratio of Total Tests to Positive Results from \code{easy_tool} Objects
#'
#' @description \code{ntpp.easy_tool} computes the ratio of the total number of
#' tests performed per positive test result from \code{easy_tool}-class objects.
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
#' The anticipated number of tests required to detect a single positive
#' \emph{nntp} is given by
#' \deqn{nntp = (SeP + (1 - Sp)(1 - P)) / SeP}
#' where \emph{Se} is sensitivity, \emph{P} is prevalence and \emph{Sp} is
#' specificity. The anticipated prevalence among those screened out is given by
#' \deqn{Puntested = ((1 - Se)P) / ((1 - Se)P + Sp (1 - P))}
#'
#' @seealso
#' \code{\link{ntpp.lasso_screenr}}
#' \code{\link{ntpp.logreg_screenr}}
#' \code{\link{ntpp.data.frame}}
#' \code{\link{ntpp.simple_screenr}}
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
#' Compute the Ratio of Total Tests to Positive Results from \code{lasso_screenr} Objects
#'
#' @description \code{ntpp.lasso_screenr} computes the ratio of the total number of
#' tests performed per positive test result from \code{lasso_screenr}-class objects.
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
#' The anticipated number of tests required to detect a single positive
#' \emph{nntp} is given by
#' \deqn{nntp = (SeP + (1 - Sp)(1 - P)) / SeP}
#' where \emph{Se} is sensitivity, \emph{P} is prevalence and \emph{Sp} is
#' specificity. The anticipated prevalence among those screened out is given by
#' \deqn{Puntested = ((1 - Se)P) / ((1 - Se)P + Sp (1 - P))}
#'
#' @examples
#' attach(uniobj1)
#' ntpp(uniobj1)
#'
#' @export
ntpp.lasso_screenr <- function(object, ..., model = c("minAIC", "minBIC"),
                               type = c("cvResults", "isResults"),
                               prev = NULL) {
     if(!class(object) == "lasso_screenr")
         stop("object not of class lasso_screenr")
     model <- match.arg(model)
     type <- match.arg(type )
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
#' Compute the Ratio of Total Tests to Positive Results from \code{logreg_screenr} Objects
#'
#' @description \code{ntpp.logreg_screenr} computes the ratio of the total number of
#' tests performed per positive test result from \code{logreg_screenr}-class objects.
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
#' The anticipated number of tests required to detect a single positive
#' \emph{nntp} is given by
#' \deqn{nntp = (SeP + (1 - Sp)(1 - P)) / SeP}
#' where \emph{Se} is sensitivity, \emph{P} is prevalence and \emph{Sp} is
#' specificity. The anticipated prevalence among those screened out is given by
#' \deqn{Puntested = ((1 - Se)P) / ((1 - Se)P + Sp (1 - P))}
#'
#' @examples
#' attach(uniobj2)
#' ntpp(uniobj2)
#'
#' @export
ntpp.logreg_screenr <- function(object, ..., type = c("cvResults", "isResults"),
                                prev = NULL) {
     if(!class(object) == "logreg_screenr")
         stop("object not of class logreg_screenr")
     type <- match.arg(type)
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


## Function ntpp.default
##
#' Compute the Ratio of Total Tests to Positive Results from a Data Frame
#'
#' @description \code{ntpp.data.frame} computes the ratio of the total number of
#' tests performed per positive test result from data frames.
#'
#' @param se  a numeric vector of sensitivities in (0,1)
#'
#' @param sp a numeric vector of sensitivities in (0,1)
#'
#' @param prev a numeric vector of prevalences of the testing condition, in (0,1)
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
#'
#' @details
#' The anticipated number of tests required to detect a single positive
#' \emph{nntp} is given by
#' \deqn{nntp = (SeP + (1 - Sp)(1 - P)) / SeP}
#' where \emph{Se} is sensitivity, \emph{P} is prevalence and \emph{Sp} is
#' specificity. The anticipated prevalence among those screened out is given by
#' \deqn{Puntested = ((1 - Se)P) / ((1 - Se)P + Sp (1 - P))}
#'
#' @importFrom dplyr between
#' @export
ntpp.default <- function(se = NULL, sp = NULL, prev = NULL, ...){
    rangecheck <- function(x, y, z, ll = 0.00001, ul = 0.9999){
        tst <- c(dplyr::between(x, ll, ul), dplyr::between(y, ll, ul), dplyr::between(z, ll, ul))
        any(!tst)
    }
    if(rangecheck(se, sp, prev)) stop("not all se, sp and prev are in (0,1)")
    object <- data.frame(sensitivity = se, specificity = sp, prev = prev)
    nnt_(object)
}


## Function ntpp.simple_screenr
##
#' Compute the Ratio of Total Tests to Positive Results from \code{simple_screenr} Objects
#'
#' @description \code{ntpp.simple_screenr} computes the ratio of the total number of
#' tests performed per positive test result from \code{simple_screenr}-class objects.
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
#' The anticipated number of tests required to detect a single positive
#' \emph{nntp} is given by
#' \deqn{nntp = (SeP + (1 - Sp)(1 - P)) / SeP}
#' where \emph{Se} is sensitivity, \emph{P} is prevalence and \emph{Sp} is
#' specificity. The anticipated prevalence among those screened out is given by
#' \deqn{Puntested = ((1 - Se)P) / ((1 - Se)P + Sp (1 - P))}
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
