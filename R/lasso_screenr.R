#################################################################################
##     R Script: lasso_screenr.R
##
##      Package: screenr
##
##  Description: Implementation of screening based on package glmpath
##
##       Author: Steve Gutreuter
##               sgutreuter@gmail.com
#################################################################################

## Function lasso_screenr
##
#' Fitting Screening Tools Using Lasso-Like Regularization of Logistic Models
#'
#' @description
#' \code{lasso_screenr} is a convenience function which combines
#' logistic regression using \emph{L}1 regularization, \emph{k}-fold
#' cross-validation, and estimation of the receiver-operating characteristic (ROC).
#' The in-sample and out-of-sample performance is estimated from the models
#' which produced the minimum AIC and minimum BIC.  Execute
#' \code{methods(class = "lasso_screenr")} to identify available methods.
#'
#' @param formula an object of class \code{stats::formula} defining the
#' testing outcome and predictor variables.
#'
#' @param data a dataframe containing the variables defined in \verb{formula}.
#' The testing outcome must be binary (0 = no/negative, 1 = yes/positive) or
#' logical (\verb{FALSE}/\verb{TRUE}).  The the predictor variables are
#' are typically binary or logical responses to questions which may be
#' predictive of the test result, but numeric variables can also be used.
#'
#' @param Nfolds the number of folds used for \emph{k}-fold cross
#' validation. Default = 10; minimum = 2, maximum = 100.
#'
#' @param L2 (logical) switch controlling penalization using the \emph{L}2 norm of
#' the parameters.  Default: \verb{TRUE}).
#'
#' @param partial_auc either a logical \verb{FALSE} or a numeric vector of the
#' form \code{c(left, right)} where left and right are numbers in the interval
#' [0, 1] specifying the endpoints for computation of the partial area under the
#' ROC curve (pAUC). The total AUC is computed if \code{partial_auc} = \verb{FALSE}.
#' Default: \code{c(0.8, 1.0)}
#'
#' @param partial_auc_focus one of \verb{"sensitivity"} or \verb{specificity},
#' specifying for which the pAUC should be computed.  \code{partial_auc.focus} is
#' ignored if \code{partial_auc} = \verb{FALSE}.  Default: \verb{"sensitivity"}.
#'
#' @param partial_auc_correct logical value indicating whether the pAUC should be
#' transformed the interval from 0.5 to 1.0. \code{partial_auc_correct} is
#' ignored if \code{partial_auc} = \verb{FALSE}. Default: \verb{TRUE}).
#'
#' @param conf_level a number between 0 and 1 specifying the confidence level
#' for confidence intervals for the (partial)AUC. Default: 0.95.
#'
#' @param boot_n number of bootstrap replications for computation of confidence
#' intervals for the (partial)AUC. Default: 4000.
#'
#' @param standardize logical; if TRUE predictors are standardized to unit
#' variance.  Default: FALSE (sensible for binary and logical predictors).
#'
#' @param seed random number generator seed for cross-validation data splitting.
#'
#' @param ... additional arguments passed to \code{\link[glmpath]{glmpath}},
#' \code{\link[pROC]{roc}}, \code{\link[pROC]{auc}} or \code{\link[pROC]{ci}} .
#'
#' @details
#' The results provide information from
#' which to choose a probability threshold above which individual out-of-sample
#' probabilies indicate the need to perform a diagnostic test.  Out-of-sample
#' performance is estimated using \emph{k}-fold cross validation.
#'
#' \code{lasso_screenr} uses the \emph{L}1 path regularizer of
#' Park and Hastie (2007), as implemented in the \code{glmpath} package.
#' Park-Hastie regularization is is similar to the conventional lasso and the
#' elastic net. It differs from the lasso with the inclusion of a very small,
#' \emph{fixed} (\verb{1e-5}) penalty on the \emph{L}2 norm of the parameter
#' vector, and differs from the elastic net in that the \emph{L}2 penalty is
#' fixed.  Like the elastic net, the Park-Hastie regularization is robust to
#' highly correlated predictors. The \emph{L}2 penalization can be turned off
#' (\code{L2 = FALSE}), in which case the regularization is similar to the
#' coventional lasso. Like all \emph{L}1 regularizers, the Park-Hastie
#' algorithm automatically "deletes" covariates by shrinking their parameter
#' estimates to 0.
#'
#' The coefficients produced by \emph{L}1 regularization are biased toward
#' zero.  Therefore one might consider refitting the model selected by
#' regularization using maximum-likelihood estimation as implemented in
#' \code{logreg_screenr}.
#'
#' The receiver-operating characteristics are computed using the \code{pROC}
#' package.
#'
#' By default, the \emph{partial} area under the ROC curve is computed from
#' that portion of the curve for which sensitivity is in the closed interval
#' [0.8, 1.0]. However, the total AUC can be obtained using the argument
#' \code{partial_auc = FALSE}.  Partial areas can be computed for either
#' ranges of sensitivity or specificity using the arguments
#' \code{partial_auc_focus} and \code{partial_auc}.  By default, partial areas
#' are standardized.
#'
#' Out-of-sample performance is estimated using \emph{k}-fold cross-validation.
#' For a gentle but Python-centric introduction to \emph{k}-fold cross-validation,
#' see
#' \url{https://machinelearningmastery.com/k-fold-cross-validation/}.
#'
#' @return
#' \code{lasso_screenr} returns (invisibly) an object of class \code{lasso_screenr}
#' containing the components:
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
#' \item{\code{RNG}}{Specification of the random-number generator used for
#' k-fold data splitting.}
#' \item{\code{RNGseed}}{RNG seed.}
#' \item{\code{cvResults}}{A list structure containing the results of \emph{k}-
#' fold cross-validation estimation of out-of-sample performance.}
#' }
#'
#' The list elements of \code{cvResutls} are:
#' \describe{
#'     \item{\code{Nfolds}}{the number folds \emph{k}}
#'     \item{\code{X_ho}}{the matrix of held-out predictors for each cross-validation
#' fold}
#'     \item{\code{minAICcvPreds}}{the held-out responses and out-of-sample predicted
#' probabilities from AIC-best model selection}
#'     \item{\code{minAICcvROC}}{the out-of-sample ROC object
#' of class \code{roc} from AIC-best model selection}
#'     \item{\code{minBICcvPreds}}{the held-out responses and out-of-sample predicted probabilities from
#' BIC-best model selection}
#'     \item{\code{minBICcvROC}}{the corresponding out-of-sample predicted probabilities
#' and ROC object from BIC-best model selection}
#' }
#'
#' @seealso \code{\link[glmpath]{glmpath}}, \code{\link[pROC]{roc}} and
#' \code{\link[pROC]{auc}}.
#'
#' @references
#' Park MY, Hastie T. \emph{L}1-regularization path algorithm for generalized linear
#' models. Journal of the Royal Statistical Society Series B. 2007;69(4):659-677.
#' \url{https://doi.org/10.1111/j.1467-9868.2007.00607.x}
#'
#' Kim J-H. Estimating classification error rate: Repeated cross-validation, repeated
#' hold-out and bootstrap. Computational Statistics and Data Analysis.
#' 2009:53(11):3735-3745. \url{http://doi.org/10.1016/j.csda.2009.04.009}
#'
#' Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez J-C,
#' Muller M. \code{pROC}: An open-source package for \code{R} and S+ to
#' analyze and compare ROC curves. BMC Bioinformatics. 2011;12(77):1-8.
#' \url{http://doi.org/10.1186/1471-2105-12-77}
#'
#' Teferi W, Gutreuter S, Bekele A et al. Adapting strategies for effective and
#' efficient pediatric HIV case finding: Risk screening tool for testing children
#' presenting at high-risk entry points. BMC Infectious Diseases. 2022; 22:480.
#' \url{http://doi.org/10.1186/s12879-022-07460-w}
#'
#' @examples
#' \dontrun{
#' data(unicorns)
#' uniobj1 <- lasso_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
#'                           data = unicorns, Nfolds = 10)
#' methods(class = class(uniobj1))
#' summary(uniobj1)
#' }
#'
#' @import pROC
#' @importFrom stats model.frame
#' @import glmpath
#' @importFrom stringr str_split
#' @export
lasso_screenr <- function(formula, data = NULL, Nfolds = 10, L2 = TRUE,
                          partial_auc = c(0.8, 1.0),
                          partial_auc_focus = "sensitivity",
                          partial_auc_correct = TRUE,
                          boot_n = 4000, conf_level = 0.95,
                          standardize = FALSE,
                          seed = Sys.time(), ... ){
    if(!inherits(formula, "formula")) stop("Specify a model formula")
    if(!is.data.frame(data)) stop("Specify a dataframe")
    if(!(Nfolds > 1 & Nfolds <= 100))
        stop("Nfolds must be in the closed interval [2, 100]")
    if(Nfolds < 5)
        warning("Nfolds < 5 is not recommended; consider this testing mode.")
    call <- match.call()
    mf <- stats::model.frame(formula, data)
    y <- as.numeric(mf[, 1])
    xclass <- lapply(mf[, -1], class)
    if(any(!(xclass %in% c("numeric", "logical", "integer"))))
        stop("The predictor variables must be numeric or logical")
    x <- apply(as.matrix(mf[, -1]), 2, as.numeric)
    N <- nrow(x)
    prev <- mean(y)
    lam2 <- ifelse(L2 == FALSE, 0, 1e-5)
    res <- glmpath::glmpath(x, y,
                            standardize = standardize,
                            family = "binomial",
                            lambda2 = lam2, ...)
    cat("\nRegularization completed.\n")
    sumry <- summary(res)
    ii <- res$new.A
    ii[length(ii)] <- TRUE
    st <- which(ii)
    pROC <- data.frame(NULL)
    message(paste0("\nEstimating in-sample (partial)AUC values along the regularization path...\n"))
    for(i in st) {
        phat <- as.vector(predict(res, newx = x, newy = y, s = i, type = "response"))
        rocx <- pROC::roc(y, phat,
                   ci = TRUE, of = "auc", conf.level = conf_level,
                   boot.n = boot_n,
                   partial.auc = partial_auc,
                   partial.auc.focus = partial_auc_focus,
                   partial.auc.correct = partial_auc_correct)
        pROC <- rbind(pROC, c(as.numeric(rocx$auc), as.numeric(rocx$ci)[c(1, 3)]))
    }
    names(pROC) <- c("pAUC", "pAUClcl", "pAUCucl")
    sumry <- cbind(data.frame(Step = st), sumry, pROC )
    rownames(sumry) <-  NULL
    AIC <- data.frame(NULL)
    BIC <- data.frame(NULL)
    pAUC <- data.frame(NULL )
    minAIC <- NULL
    minBIC <- NULL
    for(i in c("AIC", "BIC")){
        minIC <- min(sumry[[i]])
        step <- sumry[sumry[[i]] == minIC, ]$Step
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
        isROC <- pROC::roc(isPreds$y, isPreds$pred_prob, auc = TRUE,
                           ci = TRUE, of = "auc",
                           boot.n = boot_n,
                           conf.level = conf_level,
                           partial.auc = partial_auc,
                           partial.auc.focus = partial_auc_focus,
                           partial.auc.correct = partial_auc_correct)
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
    minAICcvCoef <- data.frame(NULL)
    minBICcvCoef <- data.frame(NULL)
    X_ho <- data.frame(NULL)
    message(paste0("\nPerforming ", Nfolds,"-fold cross-validation...\n"  ))
    for(j in 1:Nfolds){
        yj <- y[-holdouts[[j]]]
        xhoj <- data.frame(fold = rep(j, length(holdouts[[j]])),
                           x[holdouts[[j]],])
        rescv <- glmpath::glmpath(x[-holdouts[[j]], ], yj,
                                  standardize = standardize,
                                  family = "binomial",
                                  lambda2 = lam2, ...)
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
                             auc = TRUE, ci = TRUE, of = "auc",
                             conf.level = conf_level,
                             boot.n = boot_n,
                             partial.auc = partial_auc,
                             partial.auc.focus = partial_auc_focus,
                             partial.auc.correct = partial_auc_correct)
    minBICcvROC <- pROC::roc(minBICcvPreds$y, minBICcvPreds$pred_prob,
                             auc = TRUE, ci = TRUE, of = "auc",
                             conf.level = conf_level,
                             boot.n = boot_n,
                             partial.auc = partial_auc,
                             partial.auc.focus = partial_auc_focus,
                             partial.auc.correct = partial_auc_correct)
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
    class(result) <- "lasso_screenr"
    invisible(result)
}
