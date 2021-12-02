## Function easyTool
##
#' \code{easyTool} creates the components of an easy-to-use screening tool based
#' on the results of either \code{glmpathScreenr} or \code{logisticScreenr}.
#'
#' @param object an object of class \code{glmpathScreenr} or
#' \code{logisticScreenr}.
#'
#' @param max (numeric) the desired maximum value for the response-question
#' weights (default is 3).
#'
#' @param model (for \verb{glmpathScreenr} objects only) the desired basis
#' model. Valid options are \verb{"minAIC"} (the default) and \verb{"minBIC"}.
#'
#' @param crossval a (logical) indicator for cross-validated (\verb{TRUE}) or
#' in-sample (\verb{FALSE}) performance evaluation.
#'
#' @return
#' \code{easyTool} returns (invisibly) an object of class \code{easyTool}
#' containing:
#' \describe{
#' \item{\code{QuestionWeights}}{Weights for the screening questions obtained
#' by rescaling the non-zero-valued logistic regression coefficients to whole
#' numbers ranging from 1 to \code{max}}
#' \item{\code{Type}}{The type of test performance evaluaion ("cross-validated"
#' or "in-sample").}
#' \item{\code{Scores}}{A data frame containing the testing outcomes
#' (\code{response}) and cross-validated scores obtained as the sums of the
#' weighted responses to the set of screening questions (\code{score}).}
#' \item{\code{ROC}}{An object of class \verb{roc} containing the
#' receiver-operating characteristic produced by \code{\link[pROC]{roc}}.}
#' }
#'
#' @details
#' The estimates of the coefficients for the screening questions are rescaled to
#' whole numbers ranging from 1 to \code{max} (\code{QuestionWeights). Those are
#' used as weights for each screening question.  The cross-validation results
#' are then converted to questionnaire scores, where the score for each subject
#' is the sum of the weighted responses to each question.
#'
#' The \code{QuestionWeights} are the foundation for easy screening. For example,
#' the screening tool could consist of a simple questionnaire followed by the
#' weight for each question, expressed as a small whole number (1, ..., \code{max})
#' and/or an equal number of open circles.  The person doing the screening need
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
#' \code{coef}, \code{ntpp}, \code{plot}, \code{print} and \print{summary}
#' methods are available for \code{easyTool}-class objects.
#'
#' @seealso
#' \code{\link[screenr]{rescale_to_int}}
#'
#' \code{\link[screenr]{coef.easyTool}},\code{\link[screenr]{ntpp.easyTool}},
#' \code{\link[screenr]{plot.easyTool}}, \code{\link[screenr]{print.easyTool}}
#' and \code{\link[screenr]{summary.easyTool}}
#'
#' @examples
#' load(uniobj1)
#' tool <- easyTool(uniobj1, max = 3)
#' class(tool)
#' @export
}
easyTool <- function(object, max = 3, model = "minAIC", crossval = TRUE) {
    if(!class(object) %in% c("glmpathScreenr", "logisticScreenr"))
        stop("object not of class glmpathScreenr or logisticScreenr")
    if(class(object) == "glmpathScreenr" & !(model %in% c("minAIC", "minBIC")))
        stop("model must be one of 'minAIC' or 'minBIC'")
    if(max < 2) {
        max <- 2
        warning("max has been re-set to 2 (and 2 may not be so good)")
    }
    if(max > 10) {
        warning("max > 10; consider using the model predictions instead")
        max  <- 10
        }
    if(class(object) == "glmpathScreenr") {
        coef <- coef.glmpathScreenr(object, intercept = FALSE)
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
            X_ <- object$Smat
            parms <- coef
            resp <- object$isResults[[model]][["Preds"]][["y"]]
            fold <- NULL
        }
    } else {
        coef <- coef.logisticScreenr(object, intercept = FALSE)
        if(crossval == TRUE ){
            Nfolds <- object$Nfolds
            X_ho <- object$X_ho
            parms <- as.matrix(object$CVcoef[, -2])
            resp <- object$CVpreds$y
            fold <- object$CVpreds$fold
        } else {
            r_ <- as.character(object$formula)[2]
            d_ <- object$ModelFit$data
            X_ <- as.matrix(d_[ , !(names(x) %in% r_ )])
            parms <- coef
            resp <- d_[[r_]]
            fold <- NULL
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
        y_ <- resp
        scframe <- data.frame(response = y_, score = scores)
        attr(scframe, "Description") <- "Observed test results and scores"
        attr(scframe, "Type") <- "In-sample scores"
    }
    ROC <- pROC::roc(response ~ score, data = scframe)
    Type <- ifelse(crossval == TRUE, "cross-validated", "in-sample")
    result <- list(QuestionWeights = QuestionWeights,
                   Type = Type,
                   Scores = scframe,
                   ROC = ROC)
    attr(result, "Description")  <- "Question weights, testing results and scores from whole-numbered weighting of question responses"
    class(result) <- "easyTool"
    invisible(result)
}

## Test four cases
debugonce(easyTool)

wtf1 <- easyTool(uniobj1, max = 3)
str(wtf1)

wtf2 <- easyTool(uniobj1, max = 3, crossval = FALSE )
str(wtf2)

wtf3 <- easyTool(uniobj2, max = 3)
str(wtf3)

wtf4 <- easyTool(uniobj2, max = 3, crossval = FALSE )
str(wtf4)

rm(wtf1, wtf2, wtf3, wtf4)

## Methods:

## Function ntpp.easyTool

## Function plot.easyTool

## Function summary.easyTool

## Function coef.easyTool

## Function print.easyTool
