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
#' \code{get_what} is an S3 generic function to extract components of objects
#' produced by functions in the \code{screenr} package
#'
#' @param from an object from which to extract \code{what}.
#'
#' @param what the element to extract from \code{from}.
#'
#' @param ... additional arguments.
#'
#' @seealso \code{link[screenr]{get_what.lasso_screenr}}
#' @seealso \code{link[screenr]{get_what.logreg_screenr}}
#' @seealso \code{link[screenr]{get_what.simple_screenr}}
#' @export
get_what <- function(from, what, ...) {
    UseMethod("get_what")
}


## Function get_what.easy_tool
##
#' \code{get_what.easy_tool} is an S3 method to extract components of
#' \code{easy_tool} objects.
#'
#' @param from the \code{easy_tool}-class object from which to extract
#' the component.
#'
#' @param what the (character) name of the component to extract. Valid values are
#' \verb{"Call"}, \verb{"QuestionWeights"}, \verb{"ROCci"}, \verb{"ROC"} and
#' \verb{"Scores"}.
#'
#' @param conf_level (optional) confidence level for \code{what =} \verb{ROCci}
#'
#' @param bootreps the number of bootstrap replications for estimation of
#' confidence intervals for \code{what =} \verb{"ROCci"}.
#'
#' @param se.min minimum value of sensitivity printed for
#' \code{what =} \verb{ROCci} (default = 0.7).
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
#' The following values of \code{what} return:
#' \describe{
#' \item{\verb{"Call"}}{the function call that created \code{from}.}
#' \item{\verb{"QuestionWeights"}}{the screening question weights, which
#' are the re-scaled logistic-regression coefficients.}
#' \item{\verb{ROCci}}{a data frame containing sensitivities,
#' specificities and their confidence limits, and thresholds}
#' \item{\verb{"Scores"}}{the screening scores for each subject, which are the
#' sums of the products of the binary question responses and their
#' \code{QuestionWeights}}
#' \item{\verb{"ROC"}}{the receiver-operating characteristic for the \code{Scores}}
#' }
#'
#' @examples
#' \dontrun{
#' attach(uniobj1)
#' tool <- easy_tool(uniobj1, max = 3, crossval = TRUE)
#' ## Get and print sensitivities and specificities at feasible thresholds
#' ROCci <- get_what(from = tool, what = "ROCci")
#' print(ROCci)
#' }
#' @export
get_what.easy_tool <- function(from = NULL, what = NULL, ..., bootreps = 2000,
                             conf_level = 0.95, se.min = 0.7){
    if(!"easy_tool" %in% class(from))
        stop("Object not easy_tool class")
    if(!what %in% c("QuestionWeights", "Call", "ROC", "Scores", "ROCci"))
        stop("Invalid what argument; must be one of 'Qweights', 'Call', 'ROC', 'Scores', 'ROCci'")
    res <- from[[what]]
    if(what == "ROCci") {
        roc_  <- from[["ROC"]]
        ci_ <- pROC::ci.thresholds(roc_,
                                   boot.n = bootreps,
                                   progress = "text",
                                   conf.level = conf_level,
                                   thresholds = "local maximas")
        threshold <- attr(ci_, "thresholds")
        res <- data.frame(cbind(threshold, ci_$sensitivity,
                                ci_$specificity))
        res$threshold <- round(res$threshold + 0.01)
        names(res) <- c("Threshold_score", "se.lcl", "Sensitivity",
                        "se.ucl", "sp.lcl", "Specificity", "sp.ucl")
        row.names(res) <- 1:(dim(ci_$sensitivity)[1])
        res <- res[res$Sensitivity >= se.min, c(1, 3, 2, 4, 6, 5, 7)]
        res[is.infinite(res[,1]), 1] <- 0
    }
    invisible(res)
}


## Function get_what.lasso_screenr
##
#' \code{get_what.lasso_screenr} is an S3 method to extract components of
#' \code{lasso_screenr} objects.
#'
#' @param from the \code{lasso_screenr}-class object from which to extract
#' the component.
#'
#' @param what the (character) name of the component to extract. Valid values are
#' \verb{"glmpathObj"}, \verb{ROCci}, \verb{"cvROC"} and \verb{"isROC"}.
#'
#' @param model the (character) name of the model for which the component is
#' desired.  Valid values are \verb{"minAIC"} and \verb{"minBIC"}.
#'
#' @param conf_level (optional) confidence level for \code{what =} \verb{ROCci}
#'
#' @param bootreps the number of bootstrap replications for estimation of
#' confidence intervals for \code{what =} \verb{"ROCci"}.
#'
#' @param se.min minimum value of sensitivity printed for
#' \code{what =} \verb{ROCci} (default = 0.7).
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
#' ##   feasible thresholds
#' cvROCci <- get_what(from = uniobj1,  what = "ROCci", model = "minBIC")
#' print(cvROCci)
#' }
#'
#' @export
get_what.lasso_screenr <- function(from = NULL, what = NULL, ..., model = "minAIC",
                                   conf_level = 0.95, bootreps =  2000,
                                   se.min = 0.7){
    if(!"lasso_screenr" %in% class(from))
        stop("Object not lasso_screenr class")
    if(!what %in% c("glmpathObj", "ROCci", "cvROC", "isROC"))
        stop("Invalid what argument; must be one of 'glmpathObj', 'ROCci', 'cvROC' or 'isROC'")
    if(!model %in% c("minAIC", "minBIC"))
        stop("Invalid model argument; must be one of 'minAIC' or 'minBIC'" )
    if(what == "glmpathObj") {
        res <- from[[what]]
    } else {
        if(what == "cvROC") {
            res <- from[["cvResults"]][[model]][["ROC"]]
        } else {
            if(what == "isROC") {
                res <- from[["isResults"]][["ROC"]]
            } else {
                roc_  <- from[["cvResults"]][[model]][["ROC"]]
                ci_ <- pROC::ci.thresholds(roc_,
                                           boot.n = bootreps,
                                           progress = "text",
                                           conf.level = conf_level,
                                           thresholds = "local maximas")
                threshold <- attr(ci_, "thresholds")
                res <- data.frame(cbind(threshold, ci_$sensitivity,
                                        ci_$specificity))
                names(res) <- c("Threshold_Pr", "se.lcl", "Sensitivity",
                                "se.ucl", "sp.lcl", "Specificity", "sp.ucl")
                row.names(res) <- 1:(dim(ci_$sensitivity)[1])
                res <- res[res$Sensitivity >= se.min, c(1, 3, 2, 4, 6, 5, 7)]
                res[is.infinite(res[,1]), 1] <- 0
            }
        }
    }
    invisible(res)
}

## Function get_what.logreg_screenr
##
#' \code{get_what.logreg_screenr} is an S3 method to extract components of
#' \code{logreg_screenr} objects.
#'
#' @param from the \code{logreg_screenr}-class object from which to extract
#' the component.
#'
#' @param what the (character) name of the component to extract. Valid values are
#' \verb{"ModelFit"}, \verb{"ROCci"}, \verb{"cvROC"} and \verb{"isROC"}.
#'
#' @param conf_level (optional) confidence level for \code{what =} \verb{"ROCci"}
#'
#' @param bootreps the number of bootstrap replications for estimation of
#' confidence intervals for \code{what =} \verb{"ROCci"}.
#'
#' @param se.min minimum value of sensitivity printed for
#' \code{what =} \verb{ROCci} (default = 0.7).
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
#' by \code{\link[stats]{glm}}}.
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
#' attach(uniobj2)
#' myROCci <- get_what(from = uniobj2, what = "ROCci")
#' print(myROCci)
#' }
#'
#' @export
get_what.logreg_screenr <- function(from = NULL, what = NULL, ..., conf_level = 0.95,
                                    bootreps =  2000, se.min = 0.7) {
    if(!"logreg_screenr" %in% class(from))
        stop("from not a logreg_screenr object")
    if(!what %in% c("ModelFit", "ROCci", "cvROC", "isROC"))
        stop("Invalid what argument; must be one of 'ModelFit', 'ROCci', cvROC' or 'isROC")
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
                ci_ <- pROC::ci.thresholds(roc_,
                                           boot.n = bootreps,
                                           progress = "text",
                                           conf.level = conf_level,
                                           thresholds = "local maximas")
                threshold <- attr(ci_, "thresholds")
                res <- data.frame(cbind(threshold, ci_$sensitivity,
                                        ci_$specificity))
                names(res) <- c("Threshold_Pr", "se.lcl", "Sensitivity",
                                "se.ucl", "sp.lcl", "Specificity", "sp.ucl")
                row.names(res) <- 1:(dim(ci_$sensitivity)[1])
                res <- res[res$Sensitivity >= se.min, c(1, 3, 2, 4, 6, 5, 7)]
                res[is.infinite(res[,1]), 1] <- 0
            }
        }
    }
    invisible(res)
}


## Function get_what.simple_screenr
##
#' \code{get_what.simple_screenr} is an S3 method to extract components of
#' \code{simple_screenr} objects.
#'
#' @param from the \code{simple_screenr}-class object from which to extract
#' the component.
#'
#' @param what the (character) name of the component to extract. Valid values
#' are \verb{"ROCci"} and \verb{"isROC"}.
#'
#' @param conf_level (optional) confidence level for \code{what =} \verb{"ROCci"}
#'
#' @param bootreps the number of bootstrap replications for estimation of
#' confidence intervals for \code{what =} \verb{"ROCci"}.
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
#' }
#'
#' @examples
#' \dontrun{
#' data(unicorns)
#' too_simple <- simple_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
#'                             data = unicorns)
#' roc <- get_what(from = too_simple, what = "isROC" )
#' plot(roc)
#' }
#'
#' @export
get_what.simple_screenr <- function(from = NULL, what = NULL, ...,
                                  conf_level = 0.95, bootreps =  2000) {
    if(!"simple_screenr" %in% class(from))
        stop("from not a simple_screenr object")
    if(!what %in% c("ROCci", "isROC"))
        stop("Invalid what argument; must be ROCci or isROC")
    if(what == "isROC") what <- "ISroc"
    if(what == "ISroc"){
        res <- from[[what]]
        } else {
            ci_ <- pROC::ci.thresholds(from$ISroc, boot.n = bootreps,
                                       progress = "none",
                                       conf.level = conf_level,
                                       thresholds = "local maximas")
            threshold <- as.numeric(rownames(ci_$sensitivity)) + 0.5
            res <- data.frame(cbind(threshold, ci_$sensitivity,
                                    ci_$specificity))
        names(res) <- c("threshold", "se.low", "se.median",
                        "se.high", "sp.low", "sp.median", "sp.high")
        row.names(res)  <- 1:length(threshold)
        }
    invisible(res)
}
