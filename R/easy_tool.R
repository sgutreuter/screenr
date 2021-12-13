#################################################################################
##  R CODE COLLECTION: easy_tool.R
##
##            PACKAGE: screenr
##
##        DESCRIPTION: Implementation of easy screening
##
##         WRITTEN BY: Steve Gutreuter
##                     sgutreuter@gmail.com
##
#################################################################################

## Function easy_tool
##
#' Simplifying Screening from \code{lasso_screenr} or \code{logreg_screenr} Objects
#'
#' @description
#' \code{easy_tool} rescales the model coefficients for the screening questions
#' to whole numbers ranging from 1 to \code{max} (\code{QuestionWeights}). Those are
#' used as weights for each screening question in the simplified screening tool.
#' The test screening score is the sum of the weights, which is compared with the
#' chosen threshold.
#'
#' @param object an object of class \code{lasso_screenr} or
#' \code{logreg_screenr}.
#'
#' @param max (numeric) the desired maximum value for the response-question
#' weights (default is 3).
#'
#' @param model (for \verb{lasso_screenr} objects only) the desired basis
#' model. Valid options are \verb{"minAIC"} (the default) and \verb{"minBIC"}.
#'
#' @param crossval a (logical) indicator for cross-validated (\verb{TRUE}) or
#' in-sample (\verb{FALSE}) performance evaluation.
#'
#' @param ... additional arguments passed to \code{coef.lasso_screenr} or
#' \code{coef.logreg_screenr}
#'
#' @return
#' \code{easy_tool} returns (invisibly) an object of class \code{easy_tool}
#' containing:
#' \describe{
#' \item{\code{QuestionWeights}}{Weights for the screening questions obtained
#' by rescaling the non-zero-valued logistic regression coefficients to whole
#' numbers ranging from 1 to \code{max}.}
#' \item{\code{Type}}{The type of test performance evaluaion ("cross-validated"
#' or "in-sample").}
#' \item{\code{Scores}}{A data frame containing the testing outcomes
#' (\code{response}) and cross-validated scores obtained as the sums of the
#' weighted responses to the set of screening questions (\code{score}).}
#' \item{\code{ROC}}{An object of class \verb{roc} containing the
#' receiver-operating characteristic produced by \code{`pROC::roc`}.}
#' }
#'
#' @details
#' The \code{QuestionWeights} (see Value, below) are the foundation for easy
#' screening. For example, the screening tool could consist of a simple
#' questionnaire followed by the weight for each question, expressed as a
#' small whole number (1, ..., \code{max}) and/or an equal number of open
#' circles.  The person doing the screening need
#' only circle the numerical weight and/or fill in the circles if and only if the
#' subject provides a "yes" response to a particular question.  The person doing
#' the screening then obtains the final score for that subject by adding up the
#' circled numbers or counting the total number of filled-in circles.  Testing is
#' mandatory for consenting subjects for whom that final score equals or exceeds
#' the chosen threshold based on the receiver-operating characteristics of
#' \code{CVresults}.
#'
#' The value chosen for \code{max} involves a trade-off between the ease of
#' manual scoring and the degree to which the ROC from the re-scaling matches the
#' ROC from the model. Small values of \code{max} make manual scoring easy, and
#' sufficiently large values will match the screening performance of the model
#' fit. A value of \verb{3} may be a reasonable compromise. It is prudent to
#' compare the ROCs from a few values of \code{max} with the ROC from the model
#' and base the final choice on the trade-off between ease of manual scoring and
#' the desired combination of sensitivity and specificity.
#'
#' @note Execute \code{methods(class = "easy_tool")} to see available methods.
#'
#' @seealso
#' \code{\link{rescale_to_int}}, \code{\link{ntpp.easy_tool}},
#' \code{\link{plot.easy_tool}}, \code{\link{print.easy_tool}} and
#' \code{\link{summary.easy_tool}}
#'
#' @examples
#' attach(uniobj1)
#' tool <- easy_tool(uniobj1, max = 3)
#' class(tool)
#' @export
easy_tool <- function(object, max = 3, model = "minAIC",
                     crossval = TRUE, ...) {
    if(!class(object) %in% c("lasso_screenr", "logreg_screenr"))
        stop("object not of class lasso_screenr or logreg_screenr")
    if(class(object) == "lasso_screenr" & !(model %in% c("minAIC", "minBIC")))
        stop("model must be one of 'minAIC' or 'minBIC'")
    if(max < 2) {
        max <- 2
        warning("max has been re-set to 2 (and 2 may not be so good)")
    }
    if(max > 10) {
        warning("max > 10; consider using the model predictions instead")
        max  <- 10
    }
    call <- match.call()
    if(class(object) == "lasso_screenr") {
        coef <- coef.lasso_screenr(object, intercept = FALSE, ...)
        if(model == "minAIC") {
            coef <- coef[, 1]
        } else {
            coef <- coef[, 2]
        }
        if(crossval == TRUE) {
            Nfolds <- object$cvResults$Nfolds
            X_ho <- object[["cvResults"]][["X_ho"]]
            parms <- as.matrix(object[["cvResults"]][[model]][["Coef"]][ , -2])
            resp <- object[["cvResults"]][[model]][["Preds"]][["y"]]
            fold <- object[["cvResults"]][[model]][["Preds"]][["fold"]]
        } else {
            X_ <- object$Xmat
            parms <- coef
            resp <- object$isResults[[model]][["Preds"]][["y"]]
            fold <- rep(NA, length(resp))
        }
    } else {
        coef <- coef.logreg_screenr(object, intercept = FALSE, ...)
        if(crossval == TRUE ){
            Nfolds <- object$Nfolds
            X_ho <- object$X_ho
            parms <- as.matrix(object$CVcoef[, -2])
            resp <- object$CVpreds$y
            fold <- object$CVpreds$fold
        } else {
            r_ <- as.character(object$formula)[2]
            d_ <- object$ModelFit$data
            X_ <- as.matrix(d_[ , !(names(d_) %in% r_ )])  ## Trouble with x
            parms <- coef
            resp <- d_[[r_]]
            fold <- rep(NA, length(resp))
        }
    }
    resp <- data.frame(fold = fold, resp = resp)
    QuestionWeights <- rescale_to_int(coef, max = max)
    scframe <- data.frame(NULL)
    if(crossval == TRUE){
        for(i in 1:Nfolds) {
            pv <- parms[i, -1]
            pv[pv < 0] <- 0
            Qwts <- rescale_to_int(pv, max = max)
            X_ <- X_ho[X_ho$fold == i, -1]
            X_ <- as.matrix(X_)
            scores <- X_ %*% Qwts
            y_  <- resp$resp[resp$fold == i]
            scframe <- rbind(scframe, data.frame(fold = i,
                                                 response = y_,
                                                 score = scores))
            attr(scframe, "Description") <- "Observed test results and scores"
            attr(scframe, "Type") <- "Cross-validated scores"
        }
    } else {
        pv <- as.matrix(parms, ncol =  1)
        pv[pv < 0] <- 0
        Qwts <- rescale_to_int(pv, max = max)
        scores <- X_ %*% Qwts
        y_ <- resp[["resp"]]
        scframe <- data.frame(response = y_, score = scores)
        attr(scframe, "Description") <- "Observed test results and scores"
        attr(scframe, "Type") <- "In-sample scores"
    }
    ROC <- pROC::roc(response ~ score, data = scframe)
    Type <- ifelse(crossval == TRUE, "cross-validated", "in-sample")
    result <- list(Call = call,
                   QuestionWeights = QuestionWeights,
                   Type = Type,
                   Scores = scframe,
                   ROC = ROC)
    attr(result, "Description")  <- "Question weights, testing results and scores from whole-numbered weighting of question responses"
    class(result) <- "easy_tool"
    invisible(result)
}
