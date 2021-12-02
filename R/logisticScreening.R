#################################################################################
##       R PROGRAM: logisticScreening.R
##
##         PROJECT: R functions for HIV screening tool development
##
##      WRITTEN BY: Steve Gutreuter
##                  E-mail:  sgutreuter@gmail.com
##
#################################################################################

## Function logisticScreenr
##
#' \code{logisticScreenr} estimates  model parameters and cross-validated performance of test
#' screening based on logistic regression.
#'
#' @param formula an object of class \code{\link[stats]{formula}} defining the
#' testing
#' outcome and predictor covariates, which is passed to \code{stats::glm()}.
#'
#' @param data  the "training" sample; a data frame containing the testing
#' outcome and predictive covariates to be used for testing screening.  The
#' testing outcome must be binary (0,1) indicating negative and positive test
#' results, respectively, or logical (\verb{TRUE}/\verb{FALSE}).  The covariates
#' are typically
#' binary (0 = no, 1 = yes) responses to questions which may be predictive of
#' the test result, but any numeric or factor covariates can be used.
#'
#' @param link the character-valued name of the link function for logistic
#' regression.  Choices are \verb{"logit"} (default), \verb{"cloglog"} or
#' \verb{"probit"}.
#'
#' @param Nfolds number of folds used for \emph{k}-fold cross
#' validation (default = 10, minimum = 2, maximum = 100).
#'
#' @param ... additional arguments passsed to or from other \code{stats::glm}
#' or \code{pROC::roc}.
#'
#' @return An object of class logisticScreenr containing the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{formula}}{The formula object.}
#' \item{\code{Prevalence}}{Prevalence (proportion) of the test condition in the
#' training sample.}
#' \item{\code{ModelFit}}{An object of class \verb{glm} (See \code{\link{glm}})
#' containing the results of the model fit.}
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
#' \item{code{CVcoef}}{the estimated coefficients from cross-validation}
#' \item{code{X_ho}}{{the matrix of held-out predictors for each cross-validation fold}}
#' }
#'
#' @details
#' The results provide information from
#' which to choose a probability threshold above which individual out-of-sample
#' probabilies indicate the need to perform a diagnostic test.  Out-of-sample
#' performance is estimated using \emph{k}-fold cross validation.
#'
#' The receiver operating characteristics are computed using the \code{pROC}
#' package. See References and package documentation for additional details.
#'
#' For a gentle but python-centric introduction to \emph{k}-fold cross-validation,
#' see \link{https://machinelearningmastery.com/k-fold-cross-validation/}.
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
#' data(unicorns)
#' help(unicorns)
#' uniobj2 <- logisticScreenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6,
#'                            data = unicorns, link = "logit", Nfolds = 10)
#' summary(uniobj2)
#'
#' @alias{binomialScreenr}
#' @import pROC
#' @importFrom stats binomial predict model.response
#' @export
logisticScreenr <- function(formula,
                            data = NULL,
                            link = "logit",
                            Nfolds = 10,
                            seed = Sys.time(),
                            ...){
    if(!inherits(formula, "formula")) stop("Specify an model formula")
    if(!is.data.frame(data)) stop("Provide a data frame")
    if(!link %in% c("logit", "cloglog", "probit")) stop("Invalid link")
    call <- match.call()
    m <- match(c("formula", "data"), names(call), 0L)
    mf <- call[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mfx <- model.frame(formula, data)
    x <- as.matrix(mfx[, -1])
    dat <- eval(mf, parent.frame())
    dat <- dat[complete.cases(dat), ]
    if(Nfolds > 0.20*dim(dat)[1])
        stop("Nfolds must be < 20% of number of complete observations")
    y <- stats::model.response(dat, "numeric")
    if(!all(y %in% c(0, 1))) stop("Response variable must be binary (0, 1)")
    prev <- mean(y, na.rm = TRUE)
    lrfit <- stats::glm(formula, data = dat, family = binomial(link = link))
    parmEst <- lrfit$coeff[-1]
    if(any(parmEst < 0))
        warning("Some coefficient(s) < 0; associations should be positive.")
    is.roc <- pROC::roc(lrfit$y, lrfit$fitted.values)
    N <- nrow(dat)
    holdouts <- split(sample(1:N), 1:Nfolds)
    cv.results <- data.frame(NULL)
    cv.coef <- data.frame(NULL)
    X_ho <- data.frame(NULL)
    set.seed(seed)
    for(i in 1:Nfolds){
        xhoj <- data.frame(fold = rep(i, length(holdouts[[i]])),
                           x[holdouts[[i]],])
        X_ho <- rbind(X_ho, xhoj)
        res <- stats::glm(formula, data = dat[-holdouts[[i]], ],
                          family = binomial(link = link))
        pred.prob <- inverseLink(stats::predict(res,
                                                newdata = dat[holdouts[[i]], ]),
                                 link = link)
        y <- stats::model.response(dat[holdouts[[i]], ])
        cv.results <- rbind(cv.results,
                            data.frame(cbind(fold = rep(i, length(pred.prob)),
                                             y = y,
                                             cv.pred.prob = pred.prob)))
        coef_ <- t(as.matrix(res$coefficients))
        colnames(coef_)[1] <- "Intercept"
        coef_ <- data.frame(fold = i, coef_)
        cv.coef <- rbind(cv.coef, coef_)
    }
    attr(X_ho, "Description") <- "Hold-out predictors"
    cv.roc <- pROC::roc(cv.results$y, cv.results$cv.pred.prob,
                        auc = TRUE)
    class(cv.results) <-  c("cv.predictions", "data.frame")
    result <- list(Call = call,
                   formula = formula,
                   Prevalence = prev,
                   ModelFit = lrfit,
                   ISroc = is.roc,
                   Nfolds = Nfolds,
                   CVpreds = cv.results,
                   CVroc = cv.roc,
                   CVcoef = cv.coef,
                   X_ho = X_ho)
    class(result) <- "logisticScreenr"
    invisible(result)
}


## Function coef.logisticScreenr
##
#' \code{coef.logisticScreenr} is an S3 coef method for \code{logisticScreenr} objects.
#'
#' @param x an object of class \code{logisticScreenr}
#'
#' @param intercept (logical) retain (\code{TRUE}, default) or drop
#' (\code{FALSE}) the intercept coefficients.
#'
#' @param or return odds ratios if \verb{TRUE}; logit-scale coefficients
#' are the default.
#'
#' @details
#' Extracts the estimated coefficients from \code{logisticScreenr} objects.
#'
#' @return A numeric vector containing the estimated coefficients on the logit
#' scale.
#'
#' @examples
#' load(uniobj2)
#' class(uniobj2)
#' coef(uniobj2)
#'
#' @export
coef.logisticScreenr <- function(x, intercept =  TRUE, or = FALSE){
    stopifnot(class(x) == "logisticScreenr")
    coef_ <- x[["ModelFit"]][["coefficients"]]
    if(intercept == FALSE) coef_ <- coef_[-1]
    if(or == TRUE ) coef_ <- exp(coef_)
    coef_
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
#' load(uniobj2)
#' class(uniobj2)
#' mfit <- getWhat(from = uniobj2, what = "ModelFit")
#' mfit$coefs
#'
#' @export
getWhat.logisticScreenr <- function(from = NULL, what = NULL) {
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


## Function ntpp.logisticScreenr
##
#' \code{ntpp.logisticScreenr} is a method for computation of the anticipated
#' number of tests per positive test result
#'
#' @param object a \code{logisticScreenr}-class object produced by \code{logisticScreenr}.
#'
#' @param type (character) one of \verb{"cvResults"} (the default) or
#' \verb{"isResults"} to specify \emph{k}-fold cross-validated or in-sample
#' receiver-operating characteristics, respectively.
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
#' load(uniobj2)
#' class(uniobj2)
#' ntpp(uniobj2)
#'
#' @export
ntpp.logisticScreenr <- function(object, type = "cvResults",
                                prev = NULL) {
     if(!class(object) == "logisticScreenr")
         stop("object not of class logisticScreenr")
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


## Function plot.logisticScreenr
##
#' \code{plot.logisticScreenr} is an S3 plot method for \code{logisticScreenr} objects,
#'

#' @param x an object of class \code{logisticScreenr}.
#' @param plot_ci logical indicator for plotting point-wise confidence
#' intervals at the locally maximum subset of coordinates for
#' on sensitivity and specificity (default = \verb{TRUE}). See also
#' \code{\link[pROC]{ci.thresholds}}.
#' @param print logical indicator to return a dataframe of plot points if \verb{TRUE}
#' (default = \verb{TRUE}).
#' @param conf_level confidence level in the interval (0,1). Default is 0.95
#' producing 95\% confidence intervals
#' @param bootreps number of bootstrap replications for estimation of confidence
#' (default = 2000).
#' @param ... additional arguments passed to \code{\link[pROC]{plot.roc}} and friends.
#'
#' @return This function produces a plot as a side effect and (optionally)
#' returns a dataframe dataframe containing numerical values.
#'
#' @details
#' Plot cross-validated (out-of-sample) ROC curve with pointwise confidence
#' intevals along with the overly optimistic in-sample ROC curve.
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
#' load(uniobj2)
#' class(uniobj2)
#' plot(uniobj2)
#'
#' @importFrom graphics legend plot
#' @export
plot.logisticScreenr <- function(x, plot_ci = TRUE, print = TRUE,
                              conf_level = 0.95, bootreps = 2000, ...){
    if(!class(x) == "logisticScreenr") stop("x is not a logisticScreenr object")
    stopifnot(conf_level > 0 & conf_level < 1)
    plot(x$CVroc, print.auc = TRUE, ci = FALSE, ...)
    if(plot_ci | print){
        ciplt <- pROC::ci.thresholds(x$CVroc,
                                     boot.n = bootreps,
                                     progress = "text",
                                     conf.level = conf_level,
                                     thresholds = "local maximas")
    }
    if(print){
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
    if(print) return(citable)
}


## Function predict.logisticScreenr
##
#' \code{predict.logisticScreenr} is an S3 prediction method for objects of class \code{logisticScreenr}
#'
#' @param object an object of class \code{logisticScreenr} produced by
#' \code{\link[screenr]{logisticScreenr}}.
#'
#' @param newdata new dataframe from which predicted probabilities of positive test results are
#' desired. The dataframe must contain values of the same response variables and covariates that
#' were used to obtain \code{object}.
#'
#' @return \code{predict.logisticScreenr} returns (invisibly) a dataframe augmenting
#' \code{newdata} with the predicted probabilities of positive test results \code{phat}.
#'
#' @details This method is a convenience wrapper for \code{link[stats]{predict.glm}}.
#'
#' @examples
#' load(uniobj2)
#' class(uniobj2)
#' ## Get some new observations
#' new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
#'                         testresult = c(NA, NA), Q1 = c(0, 0), Q2 = c(0, 0),
#'                         Q3 = c(0, 0), Q4 = c(0, 0), Q5 = c(0, 1), Q6 = c(0, 1 ))
#' ## Predict the probabilities of testing positive for the new subjects
#' predict(uniobj2, newdata = new_corns)
#' @export
#'
#'
predict.logisticScreenr <- function(object = NULL, newdata = NULL, ...){
    if(!is.data.frame(newdata)) stop("Specify a dataframe")
    if(!("logisticScreenr" %in% class(object))) stop("object not a logisticScreenr object")
    form <- object$formula
    rname <- as.character(form[[2]])
    nd <- newdata
    nd[is.na(nd[[rname]]), rname] <- 0
    res <- stats::predict.glm(object$ModelFit, newdata = nd, type = "response")
    res <- data.frame(newdata, data.frame(phat = res))
    invisible(res)
}


## Function print.logisticScreenr
##
#' \code{print.logisticScreenr} is an S3 print method for \code{logisticScreenr} objects.
#'
#' @param x an object of class \code{logisticScreenr}.
#'
#' @param quote logical indicator for whether or not strings should be printed.
#'
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
#' @examples
#' load(uniobj2)
#' class(uniobj2)
#' print(uniobj2)
#'
#' @export
print.logisticScreenr <- function(x, quote = FALSE, ...){
    if(!("logisticScreenr" %in% class(x))) stop("x not logisticScreenr class")
    cat("Out-of-sample sensitivity and specificity at outcome thresholds:\n")
    df_ <- data.frame(threshold = x$CVroc$thresholds,
                      sensitivity = x$CVroc$sensitivities,
                      specificity = x$CVroc$specificities)
    print(df_)
}



## Function summary.logisticScreenr
##
#' \code{summary.logisticScreenr} is an S3 print method for \code{logisticScreenr}
#' objects.
#'
#' @param object an object of class \code{logisticScreenr} produced by function
#' \code{logisticScreenr}.
#'
#' @param diagnostics a logical value; plot model diagnostics if \verb{TRUE}.
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @return Nothing.  Summaries are printed as a side effect.
#'
#' @examples
#' load(uniobj2)
#' class(uniobj2)
#' summary(uniobj2)
#' @export
summary.logisticScreenr <- function(object, diagnostics = FALSE, ...){
    if(!("logisticScreenr" %in% class(object))) stop("object not logisticScreenr class")
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



################################   END of FILE   ################################
