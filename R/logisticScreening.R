#################################################################################
##       R PROGRAM: logisticScreening.R
##
##         PROJECT: R functions for HIV screening tool development
##
##      WRITTEN BY: Steve Gutreuter
##                  E-mail:  sgutreuter@gmail.com
##
#################################################################################

## Function logisticScreener
#' Test-Screening Tool Based on Logistic Regression
#'
#' \code{logisticScreener} estimates  model parameters and cross-validated performance of test
#' screening based on logistic regression.  The results provide information from
#' which to choose a probability threshold above which individual out-of-sample
#' probabilies indicate the need to perform a diagnostic test.  Out-of-sample
#' performance is estimated using \emph{k}-fold cross validation.
#'
#' The receiver operating characteristics are computed using the \code{pROC}
#' package. See References and package documentation for additional details.
#'
#' @param formula an object of class \code{\link[stats]{formula}} defining the
#' testing
#' outcome and predictor covariates, which is passed to \code{stats::glm()}.
#' @param data  the "training" sample; a data frame containing the testing
#' outcome and predictive covariates to be used for testing screening.  The
#' testing outcome must be binary (0,1) indicating negative and positive test
#' results, respectively, or logical (\verb{TRUE}/\verb{FALSE}).  The covariates
#' are typically
#' binary (0 = no, 1 = yes) responses to questions which may be predictive of
#' the test result, but any numeric or factor covariates can be used.
#' @param link the character-valued name of the link function for logistic
#' regression.  Choices are \verb{"logit"} (default), \verb{"cloglog"} or
#' \verb{"probit"}.
#' @param Nfolds number of folds used for \emph{k}-fold cross
#' validation (default = 40).
#' @param ... additional arguments passsed to or from other \code{stats::glm}
#' or \code{pROC::roc}.
#'
#' @return An object of class logisticScreener containing the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{ModelFit}}{An object of class "glm" (See \code{\link{glm}})
#' containing the results of the model fit.}
#' \item{\code{Prevalence}}{Prevalence (proportion) of the test condition in the
#' training sample.}
#' \item{\code{ParamEst}}{A vector containing the logistic regression parameter
#' estimates.}
#' \item{\code{ISroc}}{An object of class \code{\link[pROC]{roc}} containing
#' the "in-sample" (overly-optimistic) receiver operating characteristics,
#' and additional functions for use with this object are available in the
#' \code{pROC} package.}
#' \item{\code{CVpreds}}{An object of class \code{cv.predictions} containing
#' the data and cross-validated predicted condition \code{y}.}
#' \item{\code{CVroc}}{An object of class \code{\link[pROC]{roc}} containing
#' the \emph{k}-fold cross-validated "out-of-sample" receiver operating
#' characteristics, and additional functions for use with this object are
#' available in the \code{pROC} package.}
#' }
#'
#' @seealso \code{\link[stats]{glm}}
#'
#' @references
#' Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez J-C,
#' Müller M. \code{pROC}: An open-source package for \code{R} and S+ to
#' analyze and compare ROC curves. BMC Bioinformatics. 2011;12(77):1-8.
#' \url{http://doi.org/10.1186/1471-2105-12-77}
#'
#' @examples
#' ## Evaluate the performance of screening thresholds based on a logisitc model
#'
#' data(unicorns)
#' help(unicorns)
#' unitool <- logisticScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5,
#'                              data = unicorns, link = "logit", Nfolds = 5L)
#' summary(unitool)
#' plot(unitool)
#' \dontrun{testCounts(unitool)}
#'
#' ## Example implementation of screening based on those results:
#' ## Suppose there are new observations (excluding testing) from two previously
#' ## untested unicorns:
#'
#' new <- data.frame(ID = c('"Bernie P."', '"Alice D."'), Q1 = c(0, 0), Q2 = c(0, 0),
#'                    Q3 = c(1, 0), Q4 = c(0, 0), Q5 = c(1, 0))
#' print(new)
#'
#' ## Compute point estimates of their predicted probabilities testing positive:
#' inverseLink(as.matrix(cbind(rep(1, nrow(new)), new[, 2:6])) %*%
#'                            as.matrix(unitool$ParamEst, ncol = 1), "logit")
#' ## or, more simply,
#' predict(unitool$ModelFit, newdata = new, type = "response")
#'
#' ## If p.threshold = 0.025 is chosen as the screening threshold
#' ## (sensitivity and specificity 77\% and 69\%, respectively) then "Bernie P."
#' ## would be offered testing and "Alice D." would not.
#'
#' ## In practice, the computation of the probabilities of positive test results
#' ## among newly observed individuals might be coded outside of R using, say, a
#' ## spreadsheet.
#'
#' @import pROC
#' @importFrom stats binomial predict model.response
#' @export
logisticScreener <- function(formula,
                              data = NULL,
                              link = "logit",
                              Nfolds = 40L,
                              ...){
    if(!inherits(formula, "formula")) stop("Specify an model formula")
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
        pred.prob <- inverseLink(stats::predict(res,
                                                newdata = dat[holdouts[[i]], ]),
                                 link = link)
        y <- stats::model.response(dat[holdouts[[i]], ])
        cv.results <- rbind(cv.results,
                            data.frame(cbind(fold = rep(i, length(pred.prob)),
                                             y = y,
                                             data = dat[holdouts[[i]],],
                                             cv.pred.prob = pred.prob)))
    }
    cv.roc <- pROC::roc(cv.results$y, cv.results$cv.pred.prob,
                        auc = TRUE)
    class(cv.results) <-  c("cv.predictions", "data.frame")
    result <- list(Call = call,
                   ModelFit = lrfit,
                   Prevalence = prev,
                   ParamEst = lrfit$coeff,
                   ISroc = is.roc,
                   CVpreds = cv.results,
                   CVroc = cv.roc)
    class(result) <- "logisticScreener"
    invisible(result)
}


## Function summary.logisticScreener
#' Print Summaries of \code{logisticScreener} Objects
#'
#' @param object an object of class \code{logisticScreener} produced by function
#' \code{logisticScreening}.
#' @param diagnostics a logical value; plot model diagnostics if \verb{TRUE}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Nothing.  Summaries are printed as a side effect.
#' @export
summary.logisticScreener <- function(object, diagnostics = FALSE, ...){
    if(!("logisticScreener" %in% class(object))) stop("object not logisticScreener class")
    cat("Call:\n")
    print(object$Call)
    cat("\n\nLogistic regression model summary:\n")
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


## Function plot.logisticScreener
#' Plot ROC Curves of \code{logisticScreener} Objects
#'
#' Plot cross-validated (out-of-sample) ROC curve with pointwise confidence
#' intevals, along with the overly optimistic in-sample ROC curve.
#'
#' @param x an object of class \code{logisticScreener}.
#' @param plot_ci logical indicator for plotting point-wise confidence
#' intervals at the locally maximum subset of coordinates for
#' on sensitivity and specificity (default = \verb{TRUE}). See also
#' \code{\link[pROC]{ci.thresholds}}.
#' @param print_ci logical indicator to return a dataframe of numerical values,
#' intervals  (default = \verb{TRUE}).
#' @param conf_level confidence level in the interval (0,1). Default is 0.95
#' producing 95\% confidence intervals
#' @param bootreps number of bootstrap replications for estimation of confidence
#' (default = 2000).
#' @param ... additional arguments passed to \code{\link[pROC]{plot.roc}} and friends.
#'
#' @return This function produces a plot as a side effect and (optionally)
#' returns a dataframe dataframe containing numerical values.
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
#' @importFrom graphics legend plot
#' @export
plot.logisticScreener <- function(x, plot_ci = TRUE, print_ci = TRUE,
                              conf_level = 0.95, bootreps = 2000, ...){
    if(!class(x) == "logisticScreener") stop("x is not a logisticScreener object")
    stopifnot(conf_level > 0 & conf_level < 1)
    plot(x$CVroc, print.auc = TRUE, ci = FALSE, ...)
    if(plot_ci | print_ci){
        ciplt <- pROC::ci.thresholds(x$CVroc,
                                     boot.n = bootreps,
                                     progress = "text",
                                     conf.level = conf_level,
                                     thresholds = "local maximas")
    }
    if(print_ci){
        threshold <- attr(ciplt, "thresholds")
        citable <- data.frame(cbind(threshold, ciplt$sensitivity,
                                    ciplt$specificity))
        names(citable) <- c("threshold", "se.low", "se.median",
                            "se.high", "sp.low", "sp.median", "sp.high")
        row.names(citable) <- 1:(dim(ciplt$sensitivity)[1])
    }
    if(plot_ci) plot(ciplt)
    pROC::lines.roc(x$ISroc, lty = 3)
    legend("bottomright", legend = c("cross-validated", "in-sample"),
           lty = c(1, 3), lwd = c(2, 2))
    if(print_ci) return(citable)
}


## Function print.logisticScreener
#' Print Receiver Operating Characteristics for \code{logisticScreener} Objects
#'
#' @param x an object of class \code{logisticScreener}.
#' @param quote logical indicator for whether or not strings should be printed.
#' @param ... further arguments passed to or from other methods.
#' with surrounding quotes.
#'
#' @return Nothing. Thresholds, specificities and sensitivities are printed as a
#' side effect.
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
#' @seealso \code{\link{getROC}}
#' @export
print.logisticScreener <- function(x, quote = FALSE, ...){
    if(!("logisticScreener" %in% class(x))) stop("x not logisticScreener class")
    cat("Out-of-sample sensitivity and specificity at outcome thresholds:\n")
    df_ <- data.frame(threshold = x$CVroc$thresholds,
                      sensitivity = x$CVroc$sensitivities,
                      specificity = x$CVroc$specificities)
    print(df_)
}


## Function coef.logisticScreener
#' A method to extract the estimated coefficients from \code{logisticScreener} objects
#'
#' @param x an object of class \code{logisticScreener}
#' @param Intercept (logical) retain (\code{TRUE}, default) or drop
#' (\code{FALSE}) the intercept coefficients.
#' @param OR return odds ratios if \verb{TRUE}; logit-scale coefficients
#' are the default.
#'
#' @return A numeric vector containing the estimated coefficients on the logit
#' scale.
#' @export
coef.logisticScreener <- function(x, Intercept =  TRUE, OR = FALSE){
    stopifnot(class(x) == "logisticScreener")
    res <- x[["ModelFit"]][["coefficients"]]
    if(Intercept == FALSE) res <- res[-1]
    if(OR == TRUE ) res <- exp(res)
    res
}

### TODO: write predict.logisticScreener

################################   END of FILE   ################################
