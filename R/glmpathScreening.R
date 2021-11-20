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


## Function glmpathScreener
##
#' \code{glmpathScreener} Is a test-screening tool based on L1 regularization
#' of logistic regression based on \code{\link[glmpath]{glmpath}}.
#'
#' @param formula an object of class \code{\link[stats]{formula}} defining the
#' testing outcome and predictor covariates
#' @param data a dataframe containing the variables defined in \verb{formula}
#' @param Nfolds the number of foldes used for \emph{k}-fold cross
#' validation (default = 20, minimum = 2, maximum = 100).
#' @param criterion information criterion (IC) used to select the best model
#' (\verb{"AIC"} or \verb{"BIC"})
#' @param seed RNG seed for cross-validation data splitting
#' @param ... additional arguments passed to glmpath::glmpath or pROC::roc
#'
#' @return An object of class \code{glmpathScreener} containing the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{Prevalence}}{Prevalence of the binary response variable.}
#' \item{\code{glmpathObj}}{An object of class \code{glmpath} returned by
#' \code{glmpath::glmpath}. See \code{help(glmpath)} and
#' \code{methods(class = "glmpath")}.}
#' \item{\code{isResults}}{A list structure containing the results from the two
#' model fits which produced the minimum AIC and BIC values, respectively. The
#' results consist of \code{Coefficients} (the logit-scale parameter estimates,
#' including the intercept), \code{isPreds} (the in-sample predicted
#' probabilities) and \code{isROC} (the in-sample receiver-operating
#' characteristic (ROC) of class \code{roc}).}
#' \item{\code{cvResults}}{A list structure containing the results of \emph{k}-
#' fold cross-validation estimation of out-of-sample performance.  The list
#' elements are \code{Nfolds} (the number folds \emph{k}), \code{minAICcvPreds}
#' (the responses and corresponding out-of-sample predicted probabilities from
#' AIC-best model selection), \code{minAICcvROC} (the out-of-sample ROC object
#' of class \code{roc} from AIC-best model selection), and \code{minBICcvPreds}
#' \code{minBICcvROC} (the corresponding out-of-sample predicted probabilities
#' and ROC object from BIC-best model selection.)
#' \item{\code{RNG}}{Specification of the random-number generator used for
#' k-fold data splitting.}
#' \item{\code{RNGseed}}{RNG seed.}
#' }
#'
#' @details \code{glmpathScreener} is a convenience function which integrates
#' logistic regression using the GLM path regularizer, \emph{k}-fold
#' cross-validation and estimation of the receiver-operating characteristic.
#' The in-sample and out-of-sample performance is estimated from the models
#' which produced the minimum AIC and minimum BIC.  Execute
#' \code{methods(class = "glmpathScreener")} to identify available methods.
#'
#' @seealso \code{\link[glmpath]{glmpath}}, \code{link[pROC]{roc}}
#'
#' @references
#' Park MY, Hastie T. L1-regularization path algorithm for generalized linear
#' models. Journal of the Royal Statistical Society Series B. 2007;69(4):659-677.
#'
#' Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez J-C,
#' Müller M. \code{pROC}: An open-source package for \code{R} and S+ to
#' analyze and compare ROC curves. BMC Bioinformatics. 2011;12(77):1-8.
#' \url{http://doi.org/10.1186/1471-2105-12-77}
#'
#' @examples
#'
#' @import pROC
#' @import glmpath
#' @export
glmpathScreener <- function(formula, data = NULL, Nfolds = 20,
                            seed = Sys.time(), ... ){
    if(!inherits(formula, "formula")) stop("Specify a model formula")
    if(!is.data.frame(data)) stop("Specify a dataframe")
    if(!(Nfolds > 1 & Nfolds <= 100))
        stop("Nfolds must be in the closed interval [2, 100]")
    if(Nfolds < 5)
        warning("Nfolds < 5 is not recommended; consider this testing mode.")
    call <- match.call()
    mf <- model.frame(formula, dat)
    y <- mf[, 1]
    x <- as.matrix(mf[, -1])
    N <- nrow(x)
    prev <- mean(y)
    res <- glmpath::glmpath(x, y, standardize = FALSE, family = "binomial")
    sumry <- summary(res)
    for(i in c("AIC", "BIC")){
        minIC <- min(sumry[[i]])
        row <- rownames(sumry[sumry[[i]] == minIC, ])
        step <- as.integer(unlist(str_split(row, " "))[2])
        parmEst <- predict(res,
                           newx = x,
                           newy = y,
                           s = step,
                           type = "coefficients" )
        phat <- predict(res, x, y, s = step, type = "response")
        isPreds <- data.frame(y = y, pred_prob = as.vector(phat, mode = "numeric"))
        attr(isPreds, "Description") <- "In-sample predicted probabilities"
        isROC <- pROC::roc(isPreds$y, isPreds$pred_prob, auc = TRUE)
        assign(paste0("min", i), list(Coefficients = parmEst,
                                      Preds = isPreds,
                                      ROC = isROC))
    }
    set.seed(seed)
    holdouts <- split(sample(1:N), 1:Nfolds)
    minAICcvPreds <- data.frame(NULL)
    minBICcvPreds <- data.frame(NULL)
    minAICcvROC <- data.frame(NULL)
    minBICcvROC <- data.frame(NULL)
    for(j in 1:Nfolds){
        yj <- y[-holdouts[[j]]]
        rescv <- glmpath::glmpath(x[-holdouts[[j]], ], yj,
                                  standardize = FALSE,
                                  family = "binomial")
        sumryj <- summary(rescv)
        for(i in c("AIC", "BIC")){
            minIC <- min(sumryj[[i]])
            rowj <- rownames(sumryj[sumryj[[i]] == minIC, ])
            stepj <- as.integer(unlist(str_split(rowj, " "))[2])
            phatj <- predict(rescv,
                             newx = x[holdouts[[j]], ],
                             newy = y[holdouts[[j]]],
                             s = stepj,
                             type = "response")
            cvPreds <- data.frame(fold = rep(j, length(phatj)),
                                        y = y[holdouts[[j]]],
                                        pred_prob = as.vector(phatj,
                                                              mode = "numeric"))
            name <- paste0("min", i, "cvPreds")
            assign(name, rbind(eval(as.symbol(name)), cvPreds),
                   inherits = TRUE)
        }
    }
        attr(minAICcvPreds, "Description") <- "Out-of-sample predicted probabilities from the AIC-best model"
        attr(minBICcvPreds, "Description") <- "Out-of-sample predicted probabilities from the BIC-best model"
        minAICcvROC <- pROC::roc(minAICcvPreds$y, minAICcvPreds$pred_prob,
                                 auc = TRUE)
        minBICcvROC <- pROC::roc(minBICcvPreds$y, minBICcvPreds$pred_prob,
                                 auc = TRUE)
    result <- list(
        Call = call,
        Prevalence = prev,
        glmpathObj = res,
        Summary = sumry,
        isResults = list(minAIC = minAIC,
                         minBIC = minBIC
                         ),
        cvResults = list(Nfolds = Nfolds,
                         minAIC = list(Preds = minAICcvPreds,
                                       ROC = minAICcvROC),
                         minBIC = list(Preds = minBICcvPreds,
                                       ROC = minBICcvROC)
                         ),
        RNG = RNGkind(),
        RNGseed = seed
    )
    class(result) <- "glmpathScreener"
    invisible(result)
}


## Function coef.glmpathScreener
##
#' A method to extract the coefficients from \code{glmpathScreener} objects.
#'
#' @param object an object of class \code{glmpathScreener}.
#' @param Intercept (logical) retain (\code{TRUE}, default) or drop
#' (\code{FALSE}) the intercept coefficients.
#' @param OR return odds ratios if \verb{TRUE}; logit-scale coefficients
#' are the default
#'
#' @return a dataframe containing the estimated coefficients from the AIC-
#' and BIC-best logistic regression models.
#' @examples
#' @export
coef.glmpathScreener <- function(object, Intercept = TRUE, OR = FALSE){
    if(!("glmpathScreener" %in% class(object)))
        stop("object not glmpathScreener class")
    coef <- rbind(object$isResults$minAIC$Coefficients,
                  object$isResults$minBIC$Coefficients)
    rownames(coef) <- c("AIC-best model", "BIC-best model")
    if (Intercept == FALSE) coef <- coef
    if (OR == TRUE) coef <- exp(coef)
    data.frame(coef)
}


## Function summary.glmpathScreener
##
#' An  summary of the GLM path regularizer
#'
#' @param object a glmpathScreener object
#'
#' @return a dataframe containing the summary, including the Df, Deviance,
#' AIC and BIC for each step along the GLM path for which the active set
#' changed.
#'
#' @details This is essentially a wrapper for \code{glmpath::summary.glmpath}
#' provided for \code{glmpathScreener} objects.
#'
#' @export
summary.glmpathScreener <- function(object){
    if(!("glmpathScreener" %in% class(object)))
        stop("object not glmpathScreener class")
    res <- object$Summary
    res
}


## Function print.glmpathScreener
##
#' A print method for \code{glmpathScreener} objects
#'
#' @param object an object of class \code{glmpathScreener}
#'
#' @export
print.glmpathScreener <- function(object){
    if(!("glmpathScreener" %in% class(object)))
        stop("object not glmpathScreener class")
    cat("Function call:\n")
    print(object$Call)
    cat("\nglmpath object:\n")
    glmpath::print.glmpath(object$glmpathObj)
    cat("\nPrevalence:", object$Prevalence, "\n")
}


## Function plot.glmpathScreener
##
#' \code{plot.glmpathScreener} is a plotting method for \code{glmpathScreener}
#' objects.
#' @param x an object of class \code{glmpathScreener}.
#' @param plot_ci (logical) plot confidence intervals if \verb{TRUE}.
#' @param print_ci (logical) print a table of confidence intervals if
#' \verb{TRUE}.
#' @param model (character) select either the model which produced the
#' minimum AIC (\verb{"minAIC"}) or minimum BIC (\verb{"minBIC"}).
#' @param conf_level confidence level
#' @param bootreps the number of bootstrap replications for estimation of
#' confidence intervals.
#' @param se.min minimum value of sensitivity printed in (optional) table.
#' @param ... any additional arguments passed to \code{pROC::plot.roc} or
#' \code{pROC::lines.roc}.
#' @importFrom graphics legend plot
#' @return This function produces a plot as a side effect and (optionally)
#' returns a dataframe containing sensitivities, specificities and their
#' lower and upper confidence limits for threshold values of Pr(response = 1).
#' @details \code{plot.glmpathScreener} is an enhanced convenience wrapper for
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
#' @importFrom graphics legend plot lines
#' @import pROC
#' @export
plot.glmpathScreener <- function(x, plot_ci = TRUE, print_ci = TRUE,
                                 model = "minAIC",
                                 conf_level = 0.95, bootreps = 2000,
                                 se.min = 0.8, ...){
    if(!("glmpathScreener" %in% class(x)))
            stop("Object not glmpathScreener class")
    stopifnot(conf_level > 0 & conf_level < 1)
    if(!model %in% c("minAIC", "minBIC"))
        stop("Specify 'minAIC' or 'minBIC' for model")
    cvROC <- x$cvResults[[model]][["ROC"]]
    isROC <- x$isResults[[model]][["ROC"]]

    pROC::plot.roc(cvROC, print.auc = TRUE, ci = FALSE, ...)
    if(plot_ci | print_ci){
        ciplt <- pROC::ci.thresholds(cvROC,
                                     boot.n = bootreps,
                                     progress = "text",
                                     conf.level = conf_level,
                                     thresholds = "local maximas")
    }
    if(print_ci){
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
    if(print_ci) return(citable)
}

## Function getWhat
##
#' \code{getWhat} extracts elements from objects of class \code{glmpathScreener}.
#'
#' @param from an object of class glmpathScreener.
#' @param what (character) the element to be extracted; valid values are
#' \verb{"cvPreds"} (cross-validated predicted probabilities),
#' \verb{"cvROC"} (cross validated \code{roc}-class object),
#' \verb{"isPreds"} (predicted probabilities from the training set),
#' \verb{"isROC"} (\code{roc}-class object from the training set) and
#' \verb{"glmpathObj"} (the entire \code{glmpath}-class object).
#' @param model (character) the model from which \code{what} is desired
#' (\verb{"minAIC"} or \verb{"minBIC"}).  The value of \code{model} does not
#' matter for \code{what =} \verb{glmpathObj}, but one of the two valid values
#' must be specified (yes, that is a bit weird).
#' @return the objects specified by \code{what}.
#' @details This function is used to access the specified objects for those who
#' need to perform computations not provided by the methods for class
#' \code{glmpathScreener}.
#' @export
getWhat <- function(from,
                    what = c("cvPreds", "isPreds",
                             "cvROC", "isROC", "glmpathObj"),
                    model = "minAIC"){
    if(!("glmpathScreener" %in% class(from)))
        stop("Object not glmpathScreener class")
    if(!what %in% c("cvPreds", "isPreds", "cvROC", "isROC", "glmpathObj"))
        stop("Invalid what; valid choices are 'cvPreds', 'cvROC, 'isPreds and 'isROC'.")
    if(!what == "glmpathObj"){
        if(!model %in% c("minAIC", "minBIC"))
        stop("Specify 'minAIC' or 'minBIC' for model")
        pfx <- substring(what, 1, 2)
        typ <- substring(what, 3)
        res <- paste0(pfx, "Results")
        x <- from[[res]][[model]][[typ]]
    } else {
        x <- from$glmpathObj
    }
    invisible(x)
}
