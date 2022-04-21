################################################################################
##     R Script: get_what.R
##
##      Package: screenr
##
##  Description: S3 methods for extraction of object components
##
##       Author: Steve Gutreuter
##               sgutreuter@gmail.com
################################################################################

## Generic function get_what
##
#' S3 Methods for Extraction of Object Components
#'
#' @description
#' \code{get_what} extracts components from objects.
#'
#' @param from an object from which to extract \code{what}.
#'
#' @param what the element to extract from \code{from}.
#'
#' @param ... additional arguments.
#'
#' @seealso \code{\link{get_what.easy_tool}}, \code{\link{get_what.lasso_screenr}},
#' \code{\link{get_what.logreg_screenr}} and \code{\link{get_what.simple_screenr}}.
#' @export
get_what <- function(from, what, ...) {
    UseMethod("get_what")
}


## Function get_what.easy_tool
##
#' An S3 Method for Extraction of Components from \code{easy_tool} Objects
#'
#' @description
#' \code{get_what.easy_tool} extracts components from \code{easy_tool}-class
#' objects.
#'
#' @param from the \code{easy_tool}-class object from which to extract
#' the component.
#'
#' @param what the (character) name of the component to extract. Valid values are
#' \verb{"Call"}, \verb{"QuestionWeights"}, \verb{"ROCci"}, \verb{"ROC"} and
#' \verb{"Scores"}.
#'
#' @param conf.level (optional) confidence level for \code{what =} \verb{ROCci}
#'
#' @param bootreps the number of bootstrap replications for estimation of
#' confidence intervals for \code{what =} \verb{"ROCci"}. Default: 4000.
#'
#' @param se.min minimum value of sensitivity printed for
#' \code{what =} \verb{ROCci}.  Default: 0.8.
#'
#' @param ... optional arguments to \code{get_what} methods.
#'
#' @return The selected component is returned invisibly.
#'
#' @details
#' \code{get_what} is provided to enable easy extraction of components that
#' are not provided by the \code{plot}, \code{predict}, \code{print} or
#' \code{summary} methods.
#'
#' Valid values of \code{what} are:
#' \describe{
#' \item{\verb{"Call"}}{returns the function call that created \code{from}.}
#' \item{\verb{"QuestionWeights"}}{returns the screening question weights, which
#' are the re-scaled logistic-regression coefficients.}
#' \item{\verb{ROCci}}{returns a data frame containing sensitivities,
#' specificities and their confidence limits, and thresholds}
#' \item{\verb{"Scores"}}{returns the screening scores for each subject, which
#' are the sums of the products of the binary question responses and their
#' \code{QuestionWeights}}
#' \item{\verb{"ROC"}}{returns the receiver-operating characteristic for the
#' \code{Scores}}
#' }
#'
#' @examples
#' \dontrun{
#' attach(uniobj1)
#' tool <- easy_tool(uniobj1, max = 3, crossval = TRUE)
#' ## Get and print sensitivities and specificities at thresholds for the
#' ##   local maxima of the ROC curve
#' ROCci <- get_what(from = tool, what = "ROCci")
#' print(ROCci)
#' }
#' @export
get_what.easy_tool <- function(from = NULL, what = NULL, ..., bootreps = 4000,
                             conf.level = 0.95, se.min = 0.8){
    if(!"easy_tool" %in% class(from))
        stop("Object not easy_tool class")
    if(!what %in% c("QuestionWeights", "Call", "ROC", "Scores", "ROCci"))
        stop("Invalid what argument; must be one of 'Qweights', 'Call', 'ROC', 'Scores', 'ROCci'")
    res <- from[[what]]
    if(what == "ROCci") {
        roc_  <- from[["ROC"]]
        res <- roc_ci(roc_, bootreps = bootreps, conf.level = conf.level,
                      se.min = se.min)
        res$Threshold <- round(res$Threshold + 0.01)
    }
    invisible(res)
}


## Function get_what.lasso_screenr
##
#' An S3 Method for Extraction of Components from \code{lasso_screenr} Objects
#'
#' @description
#' \code{get_what.lasso_screenr} extracts components from \code{lasso_screenr}-class
#' objects.
#'
#' @param from the \code{lasso_screenr}-class object from which to extract
#' the component.
#'
#' @param what the character-valued name of the component to extract. Valid values are
#' \verb{"glmpathObj"}, \verb{"ROCci"}, \verb{"cvROC"} and \verb{"isROC"}.
#'
#' @param model the character-valued name of the model for which the component is
#' desired.  Valid values are \verb{"minAIC"} and \verb{"minBIC"}.
#' Default: \verb{"minAIC"}.
#'
#' @param conf.level confidence level for \code{what =} \verb{"ROCci"}. Default: 0.95.
#'
#' @param bootreps the number of bootstrap replications for estimation of
#' confidence intervals for \code{what =} \verb{"ROCci"}. Default: 4000.
#'
#' @param se.min minimum value of sensitivity printed for
#' \code{what =} \verb{ROCci}.  Default: 0.8.
#'
#' @param ... optional arguments to \code{get_what} methods.
#'
#' @return The selected component is returned invisibly.
#'
#' @details
#' \code{get_what} is provided to enable easy extraction of components that are
#' not provided by the \code{coef}, \code{plot}, \code{predict}, \code{print}
#' or \code{summary} methods.
#'
#' The following values of \code{what} return:
#' \describe{
#' \item{\verb{"glmpathObj"}}{the entire \code{glmpath}-class object produced by
#' by \code{\link[glmpath]{glmpath}}}.
#' \item{\verb{ROCci}}{a data frame containing cross-validated sensitivities,
#' specificities and their confidence limits, and thresholds}.
#' \item{\verb{"cvROC"}}{the \code{roc}-class object produced by \code{\link[pROC]{roc}}
#' containing the \emph{k}-fold cross-validated receiver-operating characteristic.}
#' \item{\verb{"isROC"}}{the \code{roc}-class object produced by \code{\link[pROC]{roc}}
#' containing the in-sample (overly optimistic) receiver-operating characteristic.}
#' }
#'
#' @examples
#' \dontrun{
#' attach(uniobj1)
#' ## Plot the coefficient paths
#' pathobj <- get_what(from = uniobj1, what = "glmpathObj", model = "minAIC")
#' plot(pathobj)
#' ## Get and print cross-validated sensitivities and specificities at
#' ##   thresholds for the local maxima of the ROC curve
#' cvROCci <- get_what(from = uniobj1,  what = "ROCci", model = "minBIC")
#' print(cvROCci)
#' }
#'
#' @export
get_what.lasso_screenr <- function(from = NULL,
                                   what = c("glmpathObj", "ROCci", "cvROC", "isROC"),
                                   ...,
                                   model = c("minAIC", "minBIC"),
                                   conf.level = 0.95, bootreps =  4000,
                                   se.min = 0.8){
    what  <- match.arg(what)
    model = match.arg(model)
    if(!"lasso_screenr" %in% class(from))
        stop("Object not lasso_screenr class")
    if(what == "glmpathObj") {
        res <- from[[what]]
    } else {
        if(what == "cvROC") {
            res <- from[["cvResults"]][[model]][["ROC"]]
        } else {
            if(what == "isROC") {
                res <- from[["isResults"]][[model]][["ROC"]]
            } else {
                roc_  <- from[["cvResults"]][[model]][["ROC"]]
                res <- roc_ci(roc_, bootreps = bootreps,
                              conf.level = conf.level, se.min = se.min)
            }
        }
    }
    invisible(res)
}


## Function get_what.logreg_screenr
##
#' An S3 Method for Extraction of Components from \code{logreg_screenr} Objects
#'
#' @description
#' \code{get_what.logreg_screenr} extracts components from \code{logreg_screenr}-class
#' objects.
#'
#' @param from the \code{logreg_screenr}-class object from which to extract
#' the component.
#'
#' @param what the (character) name of the component to extract. Valid values are
#' \verb{"ModelFit"}, \verb{"ROCci"}, \verb{"cvROC"} and \verb{"isROC"}.
#'
#' @param conf.level (optional) confidence level for \code{what =} \verb{"ROCci"}.
#' Default: 0.95.
#'
#' @param bootreps the number of bootstrap replications for estimation of
#' confidence intervals for \code{what =} \verb{"ROCci"}. Default: 4000.
#'
#' @param se.min minimum value of sensitivity printed for
#' \code{what =} \verb{ROCci}. Default: 0.8.
#'
#' @param ... optional arguments to \code{get_what} methods.
#'
#' @return The selected component is returned invisibly.
#'
#' @details
#' \code{get_what} is provided to enable easy extraction of components for those
#' who wish to perform computations that are not provided by the \code{coef},
#' \code{plot}, \code{predict}, \code{print} or \code{summary} methods.
#'
#' The following values of \code{what} return:
#' \describe{
#' \item{\verb{"ModelFit"}}{the entire \code{glm}-class object produced by
#' by \code{\link[stats]{glm}}.}
#' \item{\verb{ROCci}}{a data frame containing cross-validated sensitivities,
#' specificities and their confidence limits, and thresholds.}
#' \item{\verb{"cvROC"}}{the \code{roc}-class object produced by \code{\link[pROC]{roc}}
#' containing the \emph{k}-fold cross-validated receiver-operating characteristic.}
#' \item{\verb{"isROC"}}{the \code{roc}-class object produced by \code{\link[pROC]{roc}}
#' containing the in-sample (overly optimistic) receiver-operating characteristic.}
#' }
#'
#' @examples
#' \dontrun{
#' attach(uniobj2)
#' ## Get and print cross-validated sensitivities and specificities at
#' ##   thresholds for the local maxima of the ROC curve
#' myROCci <- get_what(from = uniobj2, what = "ROCci")
#' print(myROCci)
#' }
#'
#' @export
get_what.logreg_screenr <- function(from = NULL,
                                    what = c("ModelFit", "ROCci", "cvROC", "isROC"),
                                    ...,
                                    conf.level = 0.95,
                                    bootreps =  4000, se.min = 0.8) {
    if(!"logreg_screenr" %in% class(from))
        stop("from not a logreg_screenr object")
    what <- match.arg(what)
    if(what == "ModelFit") {
        res <- from[[what]]
    } else {
        if(what == "cvROC") {
            res <- from[["CVroc"]]
        } else {
            if(what == "isROC") {
                res <- from[["isROC"]]
            } else {
                roc_  <-  from[["CVroc"]]
                res <- roc_ci(roc_, bootreps = bootreps,
                              conf.level = conf.level, se.min = se.min)
            }
        }
    }
    invisible(res)
}


## Function get_what.simple_screenr
##
#' An S3 Method for Extraction of Components from \code{simple_screenr} Objects
#'
#' @description
#' \code{get_what.simple_screenr} extracts components from \code{simple_screenr}-class
#' objects.
#'
#' @param from the \code{simple_screenr}-class object from which to extract
#' the component.
#'
#' @param what the (character) name of the component to extract. Valid values
#' are \verb{"ROCci"} and \verb{"isROC"}.
#'
#' @param conf.level (optional) confidence level for \code{what =} \verb{"ROCci"}.
#' Default: 09.5.
#'
#' @param bootreps the number of bootstrap replications for estimation of
#' confidence intervals for \code{what =} \verb{"ROCci"}. Default: 4000.
#'
#' @param se.min minimum value of sensitivity printed for
#' \code{what =} \verb{ROCci}.  Default: 0.6.
#'
#' @param ... optional arguments to \code{get_what} methods.
#'
#' @return The selected component is returned invisibly.
#'
#' @details
#' \code{get_what} is provided to enable easy extraction of components for those
#' who wish to perform computations that are not provided by the
#' \code{plot}, \code{predict}, \code{print} or \code{summary} methods.
#'
#' The following values of \code{what} return:
#' \describe{
#' \item{\verb{"isROC"}}{the \code{roc}-class object produced by \code{\link[pROC]{roc}}
#' containing the in-sample (overly optimistic) receiver-operating characteristic.}
#' \item{\verb{"ROCci"}}{a data frame containing cross-validated sensitivities,
#' specificities and their confidence limits, and thresholds.}
#' }
#'
#' @examples
#' \dontrun{
#' data(unicorns)
#' too_simple <- simple_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
#'                             data = unicorns)
#' too_simple_roc <- get_what(from = too_simple, what = "isROC" )
#' plot(too_simple_roc)
#' }
#'
#' @export
get_what.simple_screenr <- function(from = NULL, what = c("ROCci", "isROC"), ...,
                                    conf.level = 0.95, bootreps = 4000,
                                    se.min =  0.6) {
    if(!"simple_screenr" %in% class(from))
        stop("from not a simple_screenr object")
    what <- match.arg(what)
    if(what == "isROC") what <- "ISroc"
    if(what == "ISroc"){
        res <- from[[what]]
    } else {
        res <- roc_ci(from$ISroc, bootreps = bootreps,
                      conf.level = conf.level, se.min = se.min)
    }
    res$Threshold <- round(res$Threshold + 0.01)
    invisible(res)
}
