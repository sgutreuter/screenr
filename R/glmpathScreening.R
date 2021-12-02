#################################################################################
##  R CODE COLLECTION: glmpathScreening.R
##
##            PACKAGE: screenr
##
##        DESCRIPTION: Implementation of screening based on package glmpath,
##                     including S3 methods
##
##         WRITTEN BY: Steve Gutreuter
##                     sgutreuter@gmail.com
##
#################################################################################


## Function glmpathScreenr
##
#' \code{glmpathScreenr} is a used to develop test-screening tools based on \emph{L}1
#' regularization of logistic regression based on \code{\link[glmpath]{glmpath}}.
#'
#' @param formula an object of class \code{\link[stats]{formula}} defining the
#' testing outcome and predictor covariates
#' @param data a dataframe containing the variables defined in \verb{formula}
#' @param Nfolds the number of folds used for \emph{k}-fold cross
#' validation (default = 10, minimum = 2, maximum = 100).
#' @param criterion information criterion (IC) used to select the best model
#' (\verb{"AIC"} or \verb{"BIC"})
#' @param seed RNG seed for cross-validation data splitting
#' @param ... additional arguments passed to glmpath::glmpath or pROC::roc
#'
#' @return Return (invisibly) an object of class \code{glmpathScreenr} containing
#' the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{Prevalence}}{Prevalence of the binary response variable.}
#' \item{\code{glmpathObj}}{An object of class \code{glmpath} returned by
#' \code{glmpath::glmpath}. See \code{help(glmpath)} and
#' \code{methods(class = "glmpath")}.}
#' \item{\code{Xmat}}{The matrix of predictors.}
#' \item{\code{isResults}}{A list structure containing the results from the two
#' model fits which produced the minimum AIC and BIC values, respectively. The
#' results consist of \code{Coefficients} (the logit-scale parameter estimates,
#' including the intercept), \code{isPreds} (the in-sample predicted
#' probabilities) and \code{isROC} (the in-sample receiver-operating
#' characteristic (ROC) of class \code{roc}).}
#' \item{\code{cvResults}}{A list structure containing the results of \emph{k}-
#' fold cross-validation estimation of out-of-sample performance.  The list
#' elements are:}
#' \describe{
#' \item{\code{Nfolds}}{the number folds \emph{k}}
#' \item{\code{X_ho}}{the matrix of held-out predictors for each cross-validation fold}
#' \item{\code{minAICcvPreds}}{the held-out responses and out-of-sample predicted probabilities from
#' AIC-best model selection}
#' \item{\code{minAICcvROC}}{the out-of-sample ROC object
#' of class \code{roc} from AIC-best model selection}
#' \item{\code{minBICcvPreds}}{the held-out responses and out-of-sample predicted probabilities from
#' BIC-best model selection}
#' \item{\code{minBICcvROC}}{the corresponding out-of-sample predicted probabilities
#' and ROC object from BIC-best model selection}
#' }
#' \item{\code{RNG}}{Specification of the random-number generator used for
#' k-fold data splitting.}
#' \item{\code{RNGseed}}{RNG seed.}
#' }
#'
#' @details \code{glmpathScreenr} is a convenience function which integrates
#' logistic regression using GLM-path \emph{L}1 regularization, \emph{k}-fold
#' cross-validation and estimation of the receiver-operating characteristic.
#' The in-sample and out-of-sample performance is estimated from the models
#' which produced the minimum AIC and minimum BIC.  Execute
#' \code{methods(class = "glmpathScreenr")} to identify available methods.
#'
#' The \emph{L}1 path regularizer of Park and Hastie (2007) is similar
#' to the more familiar lasso and elastic net. It differs from the
#' lasso with the inclusion of a small fixed \emph{L}2 penalty, and
#' differs from the elastic net in that the \emph{L}2 penalty is
#' fixed.  Like the elastic net, the Park-Hastie \emph{L}1 path
#' regularization is robust to highly correlated predictors.
#'
#' For a gentle but Python-centric introduction to \emph{k}-fold cross-validation,
#' see \link{https://machinelearningmastery.com/k-fold-cross-validation/}.
#'
#' @seealso \code{\link[glmpath]{glmpath}}, \code{link[pROC]{roc}}
#'
#' @references
#' Park MY, Hastie T. \emph{L}1-regularization path algorithm for generalized linear
#' models. Journal of the Royal Statistical Society Series B. 2007;69(4):659-677.
#' \url{https://doi.org/10.1111/j.1467-9868.2007.00607.x}
#'
#' Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez J-C,
#' Müller M. \code{pROC}: An open-source package for \code{R} and S+ to
#' analyze and compare ROC curves. BMC Bioinformatics. 2011;12(77):1-8.
#' \url{http://doi.org/10.1186/1471-2105-12-77}
#'
#' @examples
#' data(unicorns)
#' help(unicorns)
#' uniobj1 <- glmpathScreenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6,
#'                             data = unicorns, Nfolds = 10)
#' summary(uniobj1)
#'
#' @import pROC
#' @import glmpath
#' @export
glmpathScreenr <- function(formula, data = NULL, Nfolds = 10,
                           seed = Sys.time(), ... ){
    if(!inherits(formula, "formula")) stop("Specify a model formula")
    if(!is.data.frame(data)) stop("Specify a dataframe")
    if(!(Nfolds > 1 & Nfolds <= 100))
        stop("Nfolds must be in the closed interval [2, 100]")
    if(Nfolds < 5)
        warning("Nfolds < 5 is not recommended; consider this testing mode.")
    call <- match.call()
    mf <- model.frame(formula, data)
    y <- mf[, 1]
    x <- as.matrix(mf[, -1])
    N <- nrow(x)
    prev <- mean(y)
    res <- glmpath::glmpath(x, y, standardize = FALSE, family = "binomial")
    sumry <- summary(res)
    for(i in c("AIC", "BIC")){
        minIC <- min(sumry[[i]])
        row <- rownames(sumry[sumry[[i]] == minIC, ])
        step <- as.integer(unlist(stringr::str_split(row, " "))[2])
        parmEst <- predict(res,
                           newx = x,
                           newy = y,
                           s = step,
                           type = "coefficients")
        if(any(parmEst[-1] < 0))
            warning("Some coefficient(s) < 0; associations should be positive.")
        phat <- predict(res, x, y, s = step, type = "response")
        isPreds <- data.frame(y = y, pred_prob = as.vector(phat, mode = "numeric"))
        attr(isPreds, "Description") <- "In-sample predicted probabilities"
        isROC <- pROC::roc(isPreds$y, isPreds$pred_prob, auc = TRUE)
        assign(paste0("min", i   ), list(Coefficients = parmEst,
                                      Preds = isPreds,
                                      ROC = isROC))
    }
    set.seed(seed)
    holdouts <- split(sample(1:N), 1:Nfolds)
    minAICcvPreds <- data.frame(NULL)
    minBICcvPreds <- data.frame(NULL)
    minAICcvROC <- data.frame(NULL)
    minBICcvROC <- data.frame(NULL)
    minAICcvCoef <- data.frame(NULL)
    minBICcvCoef <- data.frame(NULL)
    X_ho <- data.frame(NULL)
    for(j in 1:Nfolds){
        yj <- y[-holdouts[[j]]]
        xhoj <- data.frame(fold = rep(j, length(holdouts[[j]])),
                           x[holdouts[[j]],])
        rescv <- glmpath::glmpath(x[-holdouts[[j]], ], yj,
                                  standardize = FALSE,
                                  family = "binomial")
        sumryj <- summary(rescv)
        for(i in c("AIC", "BIC")){
            minIC <- min(sumryj[[i]])
                      rowj <- rownames(sumryj[sumryj[[i]] == minIC, ])
            stepj <- as.integer(unlist(stringr::str_split(rowj, " "))[2])
            phatj <- predict(rescv,
                             newx = x[holdouts[[j]], ],
                             newy = y[holdouts[[j]]],
                             s = stepj,
                             type = "response")
            cvPreds <- data.frame(fold = rep(j, length(phatj)),
                                  y = y[holdouts[[j]]],
                                  pred_prob = as.vector(phatj,
                                                        mode = "numeric"))
            coef_ <- predict(rescv,
                             newx = x[holdouts[[j]], ],
                             newy = y[holdouts[[j]]],
                             s = stepj,
                             type = "coefficients")
            coef_ <- data.frame(fold = j, coef_)
            nameP <- paste0("min", i, "cvPreds")
            assign(nameP, rbind(eval(as.symbol(nameP)), cvPreds),
                   inherits = TRUE)
            nameC <- paste0("min", i, "cvCoef")
            assign(nameC, rbind(eval(as.symbol(nameC)), coef_), inherits = TRUE)
        }
        X_ho <- rbind(X_ho, xhoj)
    }
    attr(X_ho, "Description") <- "Hold-out predictors"
    attr(minAICcvPreds, "Description") <- "Out-of-sample predicted probabilities from the AIC-best model"
    attr(minBICcvPreds, "Description") <- "Out-of-sample predicted probabilities from the BIC-best model"
    attr(minAICcvCoef, "Description") <- "Out-of-sample coefficients from the AIC-best model"
    attr(minBICcvCoef, "Description") <- "Out-of-sample coefficients from the BIC-best model"
    minAICcvROC <- pROC::roc(minAICcvPreds$y, minAICcvPreds$pred_prob,
                             auc = TRUE)
    minBICcvROC <- pROC::roc(minBICcvPreds$y, minBICcvPreds$pred_prob,
                                 auc = TRUE)
    result <- list(
        Call = call,
        Prevalence = prev,
        formula = formula,
        glmpathObj = res,
        Summary = sumry,
        Xmat = x,
        isResults = list(minAIC = minAIC,
                         minBIC = minBIC
                         ),
        cvResults = list(Nfolds = Nfolds,
                         X_ho = X_ho,
                         minAIC = list(Preds = minAICcvPreds,
                                       ROC = minAICcvROC,
                                       Coef = minAICcvCoef),
                         minBIC = list(Preds = minBICcvPreds,
                                       ROC = minBICcvROC,
                                       Coef = minBICcvCoef)
                         ),
        RNG = RNGkind(),
        RNGseed = seed
    )
    class(result) <- "glmpathScreenr"
    invisible(result)
}


## Function coef.glmpathScreenr
##
#' A method to extract the estimated coefficients from \code{glmpathScreenr} objects.
#'
#' @param object an object of class \code{glmpathScreenr}.
#'
#' @param intercept (logical) retain (\code{TRUE}, default) or drop
#' (\code{FALSE}) the intercept coefficients.
#'
#' @param or return odds ratios if \verb{TRUE}; logit-scale coefficients
#' are the default.
#'
#' @return a \emph{p} x 2 matrix containing the estimated coefficients from the AIC-
#' and BIC-best logistic regression models, where \emph{p} is the number of coefficients.
#'
#' @examples
#' attach(uniobj1)
#' coef(uniobj)
#'
#' @export
coef.glmpathScreenr <- function(object, intercept = TRUE, or = FALSE){
    if(!("glmpathScreenr" %in% class(object)))
        stop("object not glmpathScreenr class")
    coef_ <- data.frame(rbind(object$isResults$minAIC$Coefficients,
                              object$isResults$minBIC$Coefficients))
    rownames(coef_) <- c("AIC-best model", "BIC-best model")
    if (intercept == FALSE) coef_ <- coef_[, -1]
    if (or == TRUE) coef_ <- exp(coef_)
    coef_ <- t(as.matrix(coef_))
    coef_
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
#' attach(uniobj1)
#' pathobj <- getWhat(from = uniobj1, what = "glmpathObj", model = "minAIC")
#' plot(pathobj)
#' cvROC <- getWhat(from = uniobj1,  what = "cvROC", model = "minBIC")
#' plot(cvROC)
#'
#' @export
getWhat.glmpathScreenr <- function(from = NULL, what = NULL, model = NULL){
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


## Function ntpp.glmpathScreenr
##
#' \code{ntpp.glmpathScreenr} is a method for computation of the anticipated
#' number of tests per positive test result
#'
#' @param object a \code{glmpathScreenr}-class object produced by \code{glmpathScreenr}.
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
#' ntpp(uniobj1)
#'
#' @export
ntpp.glmpathScreenr <- function(object, model = "minAIC", type = "cvResults",
                                prev = NULL) {
     if(!class(object) == "glmpathScreenr")
         stop("object not of class glmpathScreenr")
     if(!type %in% c("cvResults", "isResults"))
         stop("type must be 'cvResults' or 'isResults'" )
     if(is.null(prev)) prev <- object$Prevalence
     ssp <- data.frame(sensitivity = object[[type]][[model]][["ROC"]][["sensitivities"]],
                       specificity = object[[type]][[model]][["ROC"]][["specificities"]])
     ssp <- cbind(ssp, rep(prev, dim(ssp)[1]))
     names(ssp) <- c("sensitivity", "specificity", "prev")
     result <- nnt_(ssp)
     result
}


## Function plot.glmpathScreenr
##
#' \code{plot.glmpathScreenr} is a plotting method for \code{glmpathScreenr}
#' objects.
#'
#' @param x an object of class \code{glmpathScreenr}.
#' @param plot_ci (logical) plot confidence intervals if \verb{TRUE}.
#' @param print logical indicator to return a dataframe of plot points if \verb{TRUE}
#' (default = \verb{TRUE}).
#' @param model (character) select either the model which produced the
#' minimum AIC (\verb{"minAIC"}) or minimum BIC (\verb{"minBIC"}).
#' @param conf_level confidence level
#' @param bootreps the number of bootstrap replications for estimation of
#' confidence intervals.
#' @param se.min minimum value of sensitivity printed in (optional) table.
#' @param ... any additional arguments passed to \code{pROC::plot.roc} or
#' \code{pROC::lines.roc}.
#'
#' @return This function produces a plot as a side effect and (optionally)
#' returns a dataframe containing sensitivities, specificities and their
#' lower and upper confidence limits for threshold values of Pr(response = 1).
#'
#' @details \code{plot.glmpathScreenr} is an enhanced convenience wrapper for
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
#' plot(uniobj1, model = "minAIC")
#' @importFrom graphics legend plot lines
#' @import pROC
#' @export
plot.glmpathScreenr <- function(x, plot_ci = TRUE, print = TRUE,
                                 model = "minAIC",
                                 conf_level = 0.95, bootreps = 2000,
                                 se.min = 0.8, ...){
    if(!("glmpathScreenr" %in% class(x)))
            stop("Object not glmpathScreenr class")
    stopifnot(conf_level > 0 & conf_level < 1)
    if(!model %in% c("minAIC", "minBIC"))
        stop("Specify 'minAIC' or 'minBIC' for model")
    cvROC <- x$cvResults[[model]][["ROC"]]
    isROC <- x$isResults[[model]][["ROC"]]

    pROC::plot.roc(cvROC, print.auc = TRUE, ci = FALSE, ...)
    if(plot_ci | print){
        ciplt <- pROC::ci.thresholds(cvROC,
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
    pROC::lines.roc(isROC, lty = 3)
    legend("bottomright", legend = c("cross-validated", "in-sample"),
           lty = c(1, 3), lwd = c(2, 2))
    if(print) return(citable)
}


## Function predict.glmpathScreenr
##
#' \code{predict.glmpathScreenr} is a prediction method for objects of class \code{glmpathScreenr}
#'
#' @param object an object of class \code{glmpathScreenr} produced by
#' \code{\link[screenr]{glmpathScreenr}}.
#'
#' @param newdata new dataframe from which predicted probabilities of positive test results are desired.
#' The dataframe must contain values of the same response variables and covariates that were used
#' to obtain \code{obj}.
#'
#' @return \code{predict.glmpathScreenr} returns (invisibly) a dataframe augmenting the complete cases
#' in \code{newdata} with the predicted probabilities of positive test results \code{phat_minAIC} and
#' \code{phat_minBIC} from the models that produced the minimum AIC and BIC, respectively.
#'
#' @details This method is a convenience wrapper for \code{link[glmpath]{predict.glmpath}}.
#'
#' @import glmpath
#' @examples
#' attach(uniobj1)
#' ## Get some new observations
#' new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
#'                         testresult = c(NA, NA),
#'                         Q1 = c(0, 0), Q2 = c(0, 0), Q3 = c(0, 1), Q4 = c(0, 0),
#'                         Q5 = c(0, 1), Q6 = c(1, 1))
#' ## Predict the probabilities of testing positive for the new subjects
#' predict(uniobj, new_corns)
#' @export
predict.glmpathScreenr <- function(object =  NULL, newdata = NULL){
    if(!is.data.frame(newdata)) stop("Specify a dataframe")
    if(!("glmpathScreenr" %in% class(object))) stop("obj not a glmpathScreenr object")
    form  <- object$formula
    rname <- as.character(form[[2]])
    nd <- newdata
    nd[is.na(nd[[rname]]), rname] <- 0
    mf <- model.frame(form, nd)
    y <- mf[, 1]
    x <- as.matrix(mf[, -1])
    obj <- object$glmpathObj
    sAIC <- which(obj$aic == min(obj$aic))
    pAIC <- glmpath::predict.glmpath(obj, newx =  x,  newy = y,  s = sAIC, type = "response")
    sBIC <- which(obj$bic == min(obj$bic))
    pBIC <- glmpath::predict.glmpath(obj, newx =  x,  newy = y,  s = sBIC, type = "response")
    dat <- data.frame(cbind(y, x))
    names(dat)[1] <- rname
    dat[[rname]] <- newdata[[rname]]
    res <- data.frame(pAIC, pBIC)
    names(res) <- c("phat_minAIC", "phat_minBIC")
    res <- cbind(dat, res)
    invisible(res)
}


## Function print.glmpathScreenr
##
#' code{print.glmpathScreenr} is a print method for \code{glmpathScreenr} objects
#'
#' @param object an object of class \code{glmpathScreenr}
#'
#' @examples
#' attach(uniobj1)
#' print(uniobj1)
#' @export
print.glmpathScreenr <- function(object){
    if(!("glmpathScreenr" %in% class(object)))
        stop("object not glmpathScreenr class")
    cat("Function call:\n")
    print(object$Call)
    cat("\nglmpath object:\n")
    glmpath::print.glmpath(object$glmpathObj)
    cat("\nPrevalence:", object$Prevalence, "\n")
}


## Function summary.glmpathScreenr
##
#' \code{summary.glmpathScreenr} returns a summary of the GLM path regularizer
#'
#' @param object a glmpathScreenr object
#'
#' @return a dataframe containing the summary, including the Df, Deviance,
#' AIC and BIC for each step along the GLM path for which the active set
#' changed.
#'
#' @details This is essentially a wrapper for \code{glmpath::summary.glmpath}
#' provided for \code{glmpathScreenr} objects.
#' @examples
#' attach(uniobj1)
#' summary(uniobj1)
#' @export
summary.glmpathScreenr <- function(object){
    if(!("glmpathScreenr" %in% class(object)))
        stop("object not glmpathScreenr class")
    res <- object$Summary
    res
}
