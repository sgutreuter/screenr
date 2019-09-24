#################################################################################
##       R PROGRAM: binomialScreening.R
##
##         PROJECT: R fuctions for HIV screening tool development
##
##      WRITTEN BY: Steve Gutreuter, CDC/CGH/DGHT/Statistics, Estimation and
##                               Modeling Team
##                  E-mail:  sgutreuter@cdc.gov
##
##      DISCLAIMER: Although this code has been used by the Centers for Disease
##                  Control & Prevention (CDC), no warranty, expressed or
##                  implied, is made by the CDC or the U.S. Government as to the
##                  accuracy and functioning of the code and related program
##                  material nor shall the fact of distribution constitute any
##                  such warranty, and no responsibility is assumed by the CDC in
##                  connection therewith.
##
#################################################################################
#' Test-Screening Tool Based on Binomial Regression
#'
#' Estimate binomial model parameters and cross-validated performance of test
#' screening based on binomial regression.  The results provide information from
#' which to choose a probability threshold above which individual out-of-sample
#' probabilies indicate the need to perform a diagnostic test.  Out-of-sample
#' performance is estimated using \emph{k}-fold cross validation.
#'
#' @param formula an object of class \code{\link{formula}}  defining the testing
#' outcome and predictor covariates, which is passed to \code{stats::glm()}.
#' @param data  the "training" sample; a data frame containing the testing
#' outcome and predictive covariates to be used for testing screening.  The
#' testing outcome must be binary (0,1) indicating negative and positive test
#' results, respectively, or logical (TRUE/FALSE).  The covariates are typically
#' binary (0 = no, 1 = yes) responses to questions which may be predictive of
#' the test result, but any numeric or factor covariates can be used.
#' @param link the character-valued name of the link function for binomial
#' regression.  Choices are "\code{logit}" (default), "\code{cloglog}" or
#' "\code{probit}".
#' @param Nfolds an integer number of folds used for \emph{k}-fold cross
#' validation (default = 20).
#' @param ... additional arguments passsed to or from other \code{stats::glm}
#' or \code{pROC::roc}.
#'
#' @return An object of class binomscreenr containing the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{ModelFit}}{An object of class "glm" (See \code{\link{stats::glm}})}
#' containing the results of the model fit.}
#' \item{\code{Prevalence}}{Prevalence of the test condition in the training sample.}
#' \item{\code{ParamEst}}{A vector containing the binomial regression parameter estimates.}
#' \item{\code{ISroc}}{A list of class "roc" (see \code{\link{pROC::roc}})
#' containing in-sample (overly optimistic) results.}
#' \item{\code{CVpreds}}{A data frame containing \emph{k}-fold cross-validation results.}
#' \item{\code{CVroc}}{A list of class "roc" (See \code{\link{pROC::roc}})
#' containing cross-validated results.}
#' }
#'
#' @seealso \code{\link{lme4::glm}}
#'
#' @examples
#' ## Evaluate the performance of screening thresholds based on a logisitc model
#'
#' data(unicorns)
#' help(unicorns)
#' unitool <- binomialScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5,
#'                              data = unicorns, link = "logit")
#' summary(unitool)
#' plot(unitool)
#' testCounts(SensSpec = unitool)
#'
#' ## Example implementation of screening based on those results
#' ## Suppose there are new observations (excluding testing) from two previously
#' ## untested unicorns:
#'
#' new <- data.frame(ID = c('"Bernie P."', '"Alice D."'), Q1 = c(0, 0), Q2 = c(0, 0),
#'                    Q3 = c(1, 0), Q4 = c(0, 0), Q5 = c(1, 0))
#' print(new)
#' ## Compute point estimates of their predicted probabilities testing positive:
#' ## Compute point estimates of their predicted probabilities testing positive:
#' inverseLink("logit",
#'             as.matrix(cbind(c(nrow(new), nrow(new)), new[, 2:6])) %*%
#'                             as.matrix(unitool$ParmEst, ncol = 1))
#' ## or, more simply,
#' predict(unitool$ModelFit, newdata = new, type = "response")
#' ## If p.threshold = 0.025 is chosen as the screening threshold
#' ## (sensitivity and specificity 77\% and 69\%, respectively) then "Bernie P."
#' ## would be offered testing and "Alice D." would not.
#'
#' ## In practice, the computation of the probabilities of positive test results
#' ## among newly observed individuals might be coded outside of R using, say, a
#' ## spreadsheet.
#'
#' @export
binomialScreening <- function(formula,
                              data = NULL,
                              link = "logit",
                              Nfolds = 20L,
                              ...){
    if(!plyr::is.formula(formula)) stop("Specify an model formula")
    if(!is.data.frame(data)) stop("Provide a data frame")
    if(!link %in% c("logit", "cloglog", "probit")) stop("Invalid link")
    call <- match.call()
    m <- match(c("formula", "data"), names(call), 0L)
    mf <- call[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    dat <- eval(mf, parent.frame())
    dat <- dat[complete.cases(dat), ]
    if(Nfolds > 0.20*dim(dat)[1])
        stop("Nfolds must be < 20% of number of complete observations")
    y <- stats::model.response(dat, "numeric")
    if(!all(y %in% c(0, 1))) stop("Response variable must be binary (0, 1)")
    prev <- mean(y, na.rm = TRUE)
    lrfit <- stats::glm(formula, data = dat, family = binomial(link = link))
    is.roc <- pROC::roc(lrfit$y, lrfit$fitted.values)
    N <- nrow(dat)
    holdouts <- split(sample(1:N), 1:Nfolds)
    cv.results <- data.frame(NULL)
    for(i in 1:Nfolds){
        res <- stats::glm(formula, data = dat[-holdouts[[i]], ],
                          family = binomial(link = link))
        pred.prob <- inverseLink(link, predict(res, newdata = dat[holdouts[[i]], ]))
        y <- model.response(dat[holdouts[[i]], ])
        cv.results <- rbind(cv.results,
                            data.frame(cbind(fold = rep(i, length(pred.prob)),
                                             y = y,
                                             data = dat[holdouts[[i]],],
                                             cv.pred.prob = pred.prob)))
    }
    cv.roc <- pROC::roc(cv.results$y, cv.results$cv.pred.prob,
                                  auc = TRUE, ci = TRUE, of = "sp",
                                  se = seq(0, 1, 0.05), ci.type = "shape")
    class(cv.results) <-  c("cv.predictions", "data.frame")
    result <- list(Call = call,
                   ModelFit = lrfit,
                   Prevalence = prev,
                   ParamEst = lrfit$coeff,
                   ISroc = is.roc,
                   CVpreds = cv.results,
                   CVroc = cv.roc)
    class(result) <- "binomscreenr"
    invisible(result)
}



#' Print Summaries of \code{binomscreenr} Objects
#'
#' @param object an object of class \code{binomscreenr} produced by function
#' \code{binomialScreening}.
#' @param diagnostics a logical value; plot model diagnostics if \code{TRUE}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Nothing.  Summaries are printed as a side effect.
#' @export
summary.binomscreenr <- function(object, diagnostics = FALSE, ...){
    if(!("binomscreenr" %in% class(object))) stop("object not binomscreenr class")
    cat("Call:\n")
    print(object$Call)
    cat("\n\nBinomial regression model summary:\n")
    print(summary(object$ModelFit))
    if(diagnostics) plot(object$ModelFit)
    cat("\nPrevalence (In-sample prevalence of condition):\n")
    print(object$Prevalence)
    cat("\nReceiver Operating Characteristics:\n")
    is.auc <- round(as.numeric(object$ISroc$auc), digits = 4)
    cat("\nIn-sample (overly optimistic) area under the curve: ",
        is.auc, "\n", sep = "")
    cv.auc <- round(as.numeric(object$CVroc$auc), digits = 4)
    cat("Out-of-sample area under the curve: ",
        cv.auc, "\n", sep = "")
}

#' Plot ROC Curves of \code{binomscreenr} Objects
#'
#' Plot cross-validated (out-of-sample) ROC curve with pointwise 95% confidence
#' intevals on specificity (gray shaded region), along with the overly optimistic
#' in-sample ROC curve.
#'
#' @param x an object of class "binomscreenr".
#' @param ... additional arguments passed to \code{\link{pROC::plot.roc}} and friends.
#'
#' @return Nothing.  This function produces a plot as a side effect.
#'
#' @references
#' Fawcett T. An introduction to ROC analysis. Pattern Recognition Letters. 2006.
#' 27(8):861-874.
#' \url{https://www.sciencedirect.com/science/article/abs/pii/S016786550500303X?via%3Dihub}
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
#' @export
plot.binomscreenr <- function(x, ...){
    if(!class(x) == "binomscreenr") stop("x is not a binomscreenr object")
    rocobj1 <- plot(x$CVroc, print.auc = TRUE, ci = TRUE, of = "sp",
                    se = seq(0, 1, 0.01), ci.type = "shape")
    rocobj2 <- lines.roc(x$ISroc, lty = 3)
    legend("bottomright", legend = c("cross-validated", "in-sample"), lty = c(1, 3),
       lwd = c(2, 2))
}


#' Print Receiver Operating Characteristics for \code{binomscreenr} Objects
#'
#' @param x an object of class \code{binomscreenr} produced by function.
#' @param ... further arguments passed to or from other methods.
#' @param quote logical, indicating whether or not strings should be printed
#' with surrounding quotes.
#'
#' @return Nothing. Thresholds, specificities and sensitivities are printed as a
#' side effect.
#'
#' @references
#' Fawcett T. An introduction to ROC analysis. Pattern Recognition Letters. 2006.
#' 27(8):861-874.
#' \url{https://www.sciencedirect.com/science/article/abs/pii/S016786550500303X?via%3Dihub}
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
#' @seealso \code{\link{getROC}}
#' @export
print.binomscreenr <- function(x, quote = FALSE, ...){
    if(!("binomscreenr" %in% class(x))) stop("x not binomscreenr class")
    cat("Out-of-sample sensitivity and specificity at outcome thresholds:\n")
    df_ <- data.frame(threshold = x$CVroc$thresholds,
                      sensitivity = x$CVroc$sensitivities,
                      specificity = x$CVroc$specificities)
    df_
}


################################   END of FILE   ################################
