################################################################################
##     R Script: getWhat.R
##
##      Package: screenr
##
##  Description: S3 methods for extraction of object components
##
##       Author: Steve Gutreuter
##               sgutreuter@gmail.com
################################################################################


## Generic function getWhat
##
#' \code{getWhat} is an S3 generic function to extract components of objects
#' produced by functions in the \code{screenr} package
#'
#' @param from an object from which to extract \code{what}.
#'
#' @param what the element to extract from \code{from}.
#'
#' @param ... additional arguments.
#'
#' @seealso \code{link[screenr]{getWhat.glmpathScreenr}}
#' @seealso \code{link[screenr]{getWhat.logisticScreenr}}
#' @seealso \code{link[screenr]{getWhat.simpleScreenr}}
#' @export
getWhat <- function(from, what, ...) {
    UseMethod("getWhat")
}


## Function getWhat.easyTool
##
#' \code{getWhat.easyTool} is an S3 method to extract components of
#' \code{easyTool} objects.
#'
#' @param from the \code{easyTool}-class object from which to extract
#' the component.
#'
#' @param what the (character) name of the component to extract. Valid values are
#' \verb{"Call"}, \verb{"QuestionWeights"}, \verb{"ROC"} and \verb{"Scores"}.
#'
#' @param ... optional arguments to \code{getWhat} methods.
#'
#' @return The selected component is returned invisibly.
#'
#' @details
#' \code{getWhat} is provided to enable easy extraction of components for those
#' who wish to perform computations that are not provided by the
#' \code{plot}, \code{predict}, \code{print} or \code{summary} methods.
#'
#' The following values of \code{what} return:
#' \describe{
#' \item{\verb{"Call"}}{the function call that created \code{from}.}
#' \item{\verb{"QuestionWeights"}}{the screening question weights, which are the re-scale
#' logistic-regression coefficients.}
#' \item{\verb{"Scores"}}{the screening scores for each subject, which are the sums of the
#' products of the binary question responses and their \code{QuestionWeights}}
#' \item{\verb{"ROC"}}{the receiver-operating characteristic for the \code{Scores}}
#' }
#'
#' @examples
#' attach(uniobj1)
#' tool <- easyTool(uniobj1, max = 3, crossval = TRUE)
#' weights <- getWhat(from = tool, what = "QuestionWeights")
#' print(weights)
#'
#' @export
getWhat.easyTool <- function(from = NULL, what = NULL, ...){
    if(!"easyTool" %in% class(from))
        stop("Object not easyTool class")
    if(!what %in% c("QuestionWeights", "Call", "ROC", "Scores"))
        stop("Invalid what argument; must be one of 'Qweights', 'Call', 'ROC', 'Scores'")
    res <- from[[what]]
    invisible(res)
}


## Function getWhat.glmpathScreenr
##
#' \code{getWhat.glmpathScreenr} is an S3 method to extract components of
#' \code{glmpathScreenr} objects.
#'
#' @param from the \code{glmpathScreenr}-class object from which to extract
#' the component.
#'
#' @param what the (character) name of the component to extract. Valid values are
#' \verb{"glmpathObj"}, \verb{"cvROC"} and \verb{"isROC"}.
#'
#' @param model the (character) name of the model for which the component is
#' desired.  Valid values are \verb{"minAIC"} and \verb{"minBIC"}.
#' #'
#' @param ... optional arguments to \code{getWhat} methods.
#'
#' @return The selected component is returned invisibly.
#'
#' @details
#' \code{getWhat} is provided to enable easy extraction of components for those
#' who wish to perform computations that are not provided by the \code{coef},
#' \code{plot}, \code{predict}, \code{print} or \code{summary} methods.
#'
#' The following values of \code{what} return:
#' \describe{
#' \item{\verb{"glmpathObj"}}{the entire \code{glmpath}-class object produced by
#' by \code{\link[glmpath]{glmpath}}}.
#' \item{\verb{"cvROC"}}{the \code{roc}-class object produced by \code{\link[pROC]{roc}}
#' containing the \emph{k}-fold cross-validated receiver-operating characteristic.}
#' \item{\verb{"isROC"}}{the \code{roc}-class object produced by \code{\link[pROC]{roc}}
#' containing the in-sample (overly optimistic) receiver-operating characteristic.}
#' }
#'
#' @examples
#' \dontrun{
#' attach(uniobj1)
#' pathobj <- getWhat(from = uniobj1, what = "glmpathObj", model = "minAIC")
#' plot(pathobj)
#' cvROC <- getWhat(from = uniobj1,  what = "cvROC", model = "minBIC")
#' plot(cvROC)
#' }
#'
#' @export
getWhat.glmpathScreenr <- function(from = NULL, what = NULL, ..., model = "minAIC"){
    if(!"glmpathScreenr" %in% class(from))
        stop("Object not glmpathScreenr class")
    if(!what %in% c("glmpathObj", "cvROC", "isROC"))
        stop("Invalid what argument; must be one of 'glmpathObj', 'cvROC' or 'isROC'")
    if(!model %in% c("minAIC", "minBIC"))
        stop("Invalid model argument; must be one of 'minAIC' or 'minBIC'" )
    if(what == "glmpathObj") {
        res <- from[[what]]
    } else {
        if(what == "cvROC") {
            res <- from[["cvResults"]][[model]][["ROC"]]
        } else {
            if(what == "isROC") res <- from[["isResults"]][["ROC"]]
        }
    }
    invisible(res)
}

## Function getWhat.logisticScreenr
##
#' \code{getWhat.logisticScreenr} is an S3 method to extract components of
#' \code{logisticScreenr} objects.
#'
#' @param from the \code{logisticScreenr}-class object from which to extract
#' the component.
#'
#' @param what the (character) name of the component to extract. Valid values are
#' \verb{"ModelFit"}, \verb{"cvROC"} and \verb{"isROC"}.
#'
#' @param ... optional arguments to \code{getWhat} methods.
#'
#' @return The selected component is returned invisibly.
#'
#' @details
#' \code{getWhat} is provided to enable easy extraction of components for those
#' who wish to perform computations that are not provided by the \code{coef},
#' \code{plot}, \code{predict}, \code{print} or \code{summary} methods.
#'
#' The following values of \code{what} return:
#' \describe{
#' \item{\verb{"ModelFit"}}{the entire \code{glm}-class object produced by
#' by \code{\link[stats]{glm}}}.
#' \item{\verb{"cvROC"}}{the \code{roc}-class object produced by \code{\link[pROC]{roc}}
#' containing the \emph{k}-fold cross-validated receiver-operating characteristic.}
#' \item{\verb{"isROC"}}{the \code{roc}-class object produced by \code{\link[pROC]{roc}}
#' containing the in-sample (overly optimistic) receiver-operating characteristic.}
#' }
#'
#' @examples
#' \dontrun{
#' attach(uniobj2)
#' class(uniobj2)
#' mfit <- getWhat(from = uniobj2, what = "ModelFit")
#' mfit$coefs
#' }
#'
#' @export
getWhat.logisticScreenr <- function(from = NULL, what = NULL, ...) {
    if(!"logisticScreenr" %in% class(from))
        stop("from not a logisticScreenr object")
    if(!what %in% c("ModelFit", "cvROC", "isROC"))
        stop("Invalid what argument; must be one of 'ModelFit', 'cvROC' or 'isROC")
    if(what == "ModelFit") {
        res <- from[[what]]
    } else {
        if(what == "cvROC") {
            res <- from[["CVroc"]]
        } else {
            if(what == "isROC") res <- from[["isROC"]]
        }
    }
    invisible(res)
}


## Function getWhat.simpleScreenr
##
#' \code{getWhat.simpleScreenr} is an S3 method to extract components of
#' \code{simpleScreenr} objects.
#'
#' @param from the \code{simpleScreenr}-class object from which to extract
#' the component.
#'
#' @param what the (optional character) name of the component to extract. The
#' only valid value is \verb{"isROC"}.
#'
#' @param ... optional arguments to \code{getWhat} methods.
#'
#' @return The selected component is returned invisibly.
#'
#' @details
#' \code{getWhat} is provided to enable easy extraction of components for those
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
#' toosimple <- simpleScreenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6,
#'                            data = unicorns)
#' roc <- getWhat(from = toosimple, what = "isROC" )
#' plot(roc)
#' }
#'
#' @export
getWhat.simpleScreenr <- function(from = NULL, what = "isROC", ...) {
    if(!"simpleScreenr" %in% class(from))
        stop("from not a simpleScreenr object")
    if(!what %in% c("isROC"))
        stop("Invalid what argument; must be isROC")
    if(what == "isROC") what <- "ISroc"
    res <- from[[what]]
    invisible(res)
}
