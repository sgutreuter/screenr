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
#' Fitting Screening Tools Using Lasso-Like Regularization of Logistic Regression
#'
#' @description
#' \code{lasso_screenr} is a convenience function which integrates
#' logistic regression using \emph{L}1 regularization, \emph{k}-fold
#' cross-validation and estimation of the receiver-operating characteristic.
#' The in-sample and out-of-sample performance is estimated from the models
#' which produced the minimum AIC and minimum BIC.  Execute
#' \code{methods(class = "lasso_screenr")} to identify available methods.
#'
#' @param formula an object of class \code{stats::formula} defining the
#' testing outcome and predictor covariates.
#'
#' @param data a dataframe containing the variables defined in \verb{formula}.
#'
#' @param Nfolds the number of folds used for \emph{k}-fold cross
#' validation (default = 10, minimum = 2, maximum = 100).
#'
#' @param L2 (logical) switch controlling penalization using the \emph{L}2 norm of the
#' parameters (default = \verb{TRUE}).
#'
#' @param seed random number generator seed for cross-validation data splitting.
#'
#' @param ... additional arguments passed to \code{\link[glmpath]{glmpath}} or
#' \code{\link[pROC]{roc}}.
#'
#' @details
#' \code{lasso_screenr} uses the \emph{L}1 path regularizer of
#' Park and Hastie (2007), which is similar
#' to conventional lasso and elastic net. It differs from the
#' lasso with the inclusion of a small fixed (\verb{1e-5}) penalty on the
#' \emph{L}2 norm of the parameters, and
#' differs from the elastic net in that the \emph{L}2 penalty is
#' fixed.  Like the elastic net, the Park-Hastie regularization is robust to
#' highly correlated predictors.
#'
#' The \emph{L}2 penalization can be turned off (\code{L2 = FALSE}), in which case
#' the regularization is similar to the coventional lasso.
#'
#' For a gentle but Python-centric introduction to \emph{k}-fold cross-validation,
#' see \url{https://machinelearningmastery.com/k-fold-cross-validation/}.
#'
#' @return
#' Return (invisibly) an object of class \code{lasso_screenr} containing
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
#' @seealso \code{\link[glmpath]{glmpath}}, \code{\link[pROC]{roc}}
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
#' uniobj1 <- lasso_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
#'                           data = unicorns, Nfolds = 10)
#' summary(uniobj1)
#'
#' @import pROC
#' @importFrom stats model.frame
#' @import glmpath
#' @importFrom stringr str_split
#' @export
lasso_screenr <- function(formula, data = NULL, Nfolds = 10, L2 = TRUE,
                           seed = Sys.time(), ... ){
    if(!inherits(formula, "formula")) stop("Specify a model formula")
    if(!is.data.frame(data)) stop("Specify a dataframe")
    if(!(Nfolds > 1 & Nfolds <= 100))
        stop("Nfolds must be in the closed interval [2, 100]")
    if(Nfolds < 5)
        warning("Nfolds < 5 is not recommended; consider this testing mode.")
    call <- match.call()
    mf <- stats::model.frame(formula, data)
    y <- mf[, 1]
    x <- as.matrix(mf[, -1])
    N <- nrow(x)
    prev <- mean(y)
    lam2 <- ifelse(L2 == TRUE, 0, 1e-5)
    res <- glmpath::glmpath(x, y,
                            standardize = FALSE,
                            family = "binomial",
                            lambda2 = lam2, ...)
    sumry <- summary(res)
    minAIC <- data.frame(NULL)
    minBIC <- data.frame(NULL)
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
    for(j in 1:Nfolds){
        yj <- y[-holdouts[[j]]]
        xhoj <- data.frame(fold = rep(j, length(holdouts[[j]])),
                           x[holdouts[[j]],])
        rescv <- glmpath::glmpath(x[-holdouts[[j]], ], yj,
                                  standardize = FALSE,
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
    class(result) <- "lasso_screenr"
    invisible(result)
}