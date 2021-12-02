#################################################################################
##  R CODE COLLECTION: easyTool.R
##
##            PACKAGE: screenr
##
##        DESCRIPTION: Implementation of easy screening
##
##         WRITTEN BY: Steve Gutreuter
##                     sgutreuter@gmail.com
##
#################################################################################


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
#' attach(uniobj1)
#' tool <- easyTool(uniobj1, max = 3)
#' class(tool)
#' @export
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
    call <- match.call()
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
            X_ <- object$Xmat
            parms <- coef
            resp <- object$isResults[[model]][["Preds"]][["y"]]
            fold <- rep(NA, length(resp))
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
    class(result) <- "easyTool"
    invisible(result)
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
#' tool <- easyTool(uniobj1, max = c, crossval = TRUE)
#' weights <- getWhat(from = tool, what = "QuestionWeights")
#' print(weights)
#'
#' @export
getWhat.easyTool <- function(from = NULL, what = NULL, model = NULL){
    if(!"easyTool" %in% class(from))
        stop("Object not easyTool class")
    if(!what %in% c("QuestionWeights", "Call", "ROC", "Scores"))
        stop("Invalid what argument; must be one of 'Qweights', 'Call', 'ROC', 'Scores'")
    res <- from[[what]]
    invisible(res)
}


## Function ntpp.easyTool
##
#' \code{ntpp.easyTool} is a method for computation of the anticipated
#' number of tests per positive test result
#'
#' @param object an \code{easyTool}-class object produced by \code{easyTool}.
#'
#' @param prev an optional prevalence proportion for the test outcome; if missing
#' the prevalence is obtained from \code{object}.
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
#' tool <- easyTool(uniobj1, max = c, crossval = TRUE)
#' ntpp(tool)
#'
#' @export
ntpp.easyTool <- function(object, prev = NULL) {
     if(!class(object) == "easyTool")
         stop("object not of class easyTool")
     if(is.null(prev)) prev <- mean(object$Scores$response, na.rm = TRUE)
     ssp <- data.frame(sensitivity = object$ROC$sensitivities,
                       specificity = object$ROC$specificities)
     ssp <- cbind(ssp, rep(prev, dim(ssp)[1]))
     names(ssp) <- c("sensitivity", "specificity", "prev")
     result <- nnt_(ssp)
     result
}


## Function plot.easyTool
##
#' \code{plot.easyTool} is a plotting method for \code{easyTool}
#' objects.
#' @param x an object of class \code{easyTool}.
#' @param plot_ci (logical) plot confidence intervals if \verb{TRUE}.
#' @param print logical indicator to return a dataframe of plot points if \verb{TRUE}
#' (default = \verb{TRUE}).
#' @param conf_level confidence level
#' @param bootreps the number of bootstrap replications for estimation of
#' confidence intervals.
#' @param se.min minimum value of sensitivity printed in (optional) table.
#' @param ... any additional arguments passed to \code{pROC::plot.roc} or
#' \code{pROC::lines.roc}.
#'
#' @importFrom graphics legend plot
#'
#' @return This function produces a plot as a side effect and (optionally)
#' returns a dataframe containing sensitivities, specificities and their
#' lower and upper confidence limits for threshold values of Pr(response = 1).
#'
#' @details \code{plot.easyTool} is an enhanced convenience wrapper for
#' \link{\code{pROC::plot.roc}}.  The table is useful for identifying the
#' minimum predicted response probabilities associated with particular
#' sensitivities.  The sensitivities and specificities are the coordinates at
#' change points in the cross-validated ROC curve, and the values of Threshold_Pr
#' are the values of lower bound of the predicted probability that acheives those
#' sensitivities and specificities.  For example, if Threshold_P = 0.002, then
#' classifiction as positive of all those subjects for whom the predicted response
#' probability is greater than or equal to 0.002 will achieve the corresponding
#' sensitivity and specificity at the prescribed confidence level in new
#' data drawn from the same distribution as the training sample.
#'
#' @references
#' Fawcett T. An introduction to ROC analysis. Pattern Recognition Letters. 2006.
#' 27(8):861-874.
#' \url{https://doi.org/10.1016/j.patrec.2005.10.010}
#'
#' Linden A. Measuring diagnostic and predictive accuracy in disease
#' management: an introduction to receiver operating characteristic (ROC) analysis.
#' Journal of Evaluation in Clinical Practice. 2006; 12(2):132-139.
#' \url{https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1365-2753.2005.00598.x}
#'
#' Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez J-C, Muller M.
#' pROC: an open-source package for R and S+ to analyze and compare ROC curves.
#' BMC Bioinformatics 2011; 12:77. \url{https://www.biomedcentral.com/1471-2105/12/77}
#'
#' @examples
#' attach(uniobj1)
#' tool <- easyTool(uniobj1, max = c, crossval = TRUE)
#' plot(tool)
#' @importFrom graphics legend plot lines
#' @import pROC
#' @export
plot.easyTool <- function(x, plot_ci = TRUE, print = TRUE,
                                 conf_level = 0.95, bootreps = 2000,
                                 se.min = 0.8, ...){
    if(!("easyTool" %in% class(x)))
            stop("Object not easyTool class")
    stopifnot(conf_level > 0 & conf_level < 1)
    roc_  <- x$ROC
    pROC::plot.roc(roc_, print.auc = TRUE, ci = FALSE, ...)
    if(plot_ci | print){
        ciplt <- pROC::ci.thresholds(roc_,
                                     boot.n = bootreps,
                                     progress = "text",
                                     conf.level = conf_level,
                                     thresholds = "local maximas")
    }
    if(print){
        threshold <- attr(ciplt, "thresholds")
        citable <- data.frame(cbind(threshold, ciplt$sensitivity,
                                    ciplt$specificity))
        names(citable) <- c("Threshold_Pr", "se.lcl", "Sensitivity",
                            "se.ucl", "sp.lcl", "Specificity", "sp.ucl")
        row.names(citable) <- 1:(dim(ciplt$sensitivity)[1])
        citable <- citable[citable$Sensitivity >= se.min, c(1, 3, 2, 4, 6, 5, 7)]
        citable[is.infinite(citable[,1]), 1] <- 0
    }
    if(plot_ci) plot(ciplt)
    if(print) return(citable)
}


## Function print.easyTool
##
#' \code{print.easyTool} is a print method for \code{easyTool} objects
#'
#' @param object an object of class \code{easyTool}
#'
#' @examples
#' load(uniobj1)
#' class(uniobj1)
#' print(uniobj1)
#' @export
print.easyTool <- function(object){
    if(!("easyTool" %in% class(object)))
        stop("object not easyTool class")
    cat("Function call:\n")
    print(object$Call)
    cat("\nType:\n")
    print(object$Type )
    cat("\nQuestion weights:\n")
    print(object$QuestionWeights )
}


## Function summary.easyTool
##
#' \code{summary.easyTool} returns a summary of the receiver-operating
#' characteristic produced by \code{easyTool}
#'
#' @param object an \code{easyTool} object
#'
#' @return a dataframe containing the summary, including the Df, Deviance,
#' AIC and BIC for each step along the GLM path for which the active set
#' changed.
#'
#' @details This is essentially a wrapper for \code{glmpath::summary.glmpath}
#' provided for \code{glmpathScreenr} objects.
#' @examples
#' load(uniobj1)
#' class(uniobj1)
#' summary(uniobj1)
#' @import pROC
#' @export
summary.easyTool <- function(object){
    if(!("easyTool" %in% class(object)))
        stop("object not easyTool class")
    print(object$ROC)
}
