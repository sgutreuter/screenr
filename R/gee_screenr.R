#################################################################################
##       R PROGRAM: gee_screenr.R
##
##         PROJECT: R function for HIV screening based on ordinary logistic
##                  regression.
##
##      WRITTEN BY: Steve Gutreuter
##                  E-mail:  sgutreuter@gmail.com
##
#################################################################################

## Function gee_screenr
##
#' Fitting Screening Tools Using GEE Estimation of Logistic Models
#'
#' @description
#' \code{gee_screenr} is a convenience function which integrates GEE estimation
#' of logsitic models, \emph{k}-fold  cross-validation and estimation of the
#' receiver-operating characteristic. GEE estimation accommodates cluster sampling.
#'
#' @param formula an object of class \code{stats::formula} defining the
#' testing outcome and predictor covariates, which is passed to \code{stats::glm()}.
#'
#' @param id a vector identifying the sampling clusters.
#'
#' @param data a dataframe containing the variables defined in \verb{formula}. The
#' testing outcome must be binary (0,1) indicating negative and positive test
#' results, respectively, or logical (\verb{TRUE}/\verb{FALSE}).  The covariates
#' are typically binary (0 = no, 1 = yes) responses to questions which may be
#' predictive of the test result, but any numeric or factor covariates can be used.
#'
#' @param link the character-valued name of the link function for logistic
#' regression.  Choices are \verb{"logit"}, \verb{"cloglog"} or
#' \verb{"probit"}. Default: \verb{"logit"}.
#'
#' @param corstr a character string specifying the correlation structure. The
#' following are permitted: \verb{"independence"}, \verb{"exchangeable"} and
#' \verb{"unstructured"}.  Default: \verb{independence}
#'
#' @param Nfolds number of folds used for \emph{k}-fold cross
#' validation (minimum = 2, maximum = 100).  Default: 10.
#'
#' @param partial_auc either a logical \verb{FALSE} or a numeric vector of the
#' form \code{c(left, right)} where left and right are numbers in the interval
#' [0, 1] specifying the endpoints for computation of the partial area under the
#' ROC curve (pAUC). The total AUC is computed if \code{partial_auc} = \verb{FALSE}.
#' Default: \code{c(0.8, 1.0)}.
#'
#' @param partial_auc_focus one of \verb{"sensitivity"} or \verb{specificity},
#' specifying for which the pAUC should be computed.  \code{partial_auc_focus} is
#' ignored if \code{partial_auc} = \verb{FALSE}.  Default: \verb{"sensitivity"}.
#'
#' @param partial_auc_correct logical value indicating whether the pAUC should be
#' transformed the interval from 0.5 to 1.0. \code{partial_auc_correct} is
#' ignored if \code{partial_auc} = \verb{FALSE}. Default: \verb{TRUE}).
#'
#' @param conf_level a number between 0 and 1 specifying the confidence level
#' for confidence intervals for the (partial)AUC. Default: 0.95.
#'
#' @param boot_n Number of bootstrap replications for computation of confidence
#' intervals for the (partial)AUC. Default: 4000.
#'
#' @param seed random-number generator seed for cross-validation data splitting.
#'
#' @param ... additional arguments passsed to or from other \code{geepack::geeglm}
#' or \code{pROC::roc}.
#'
#' @return \code{gee_screenr} returns an object of class \code{gee_screenr},
#' which inherits from class \code{logreg_screenr}, containing the elements:
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
#' \item{\code{CVcoef}}{the estimated coefficients from cross-validation}
#' \item{\code{X_ho}}{{the matrix of held-out predictors for each cross-validation fold}}
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
#' By default, the \emph{partial} area under the ROC curve is computed from
#' that portion of the curve for which sensitivity is in the closed interval
#' [0.8, 1.0]. However, the total AUC can be obtained using the argument
#' \code{partial_auc = FALSE}.  Partial areas can be computed for either
#' ranges of sensitivity or specificity using the arguments
#' \code{partial_auc_focus} and \code{partial_auc}.  By default, partial areas
#' are standardized.
#'
#' Out-of-sample performance is estimated using \emph{k}-fold cross-validation.
#' For a gentle but python-centric introduction to \emph{k}-fold cross-validation,
#' see \url{https://machinelearningmastery.com/k-fold-cross-validation/}.
#'
#' @seealso \code{\link[geepack]{geeglm}}, \code{\link[pROC]{roc}} and
#' \code{\link[pROC]{auc}}
#'
#' @references
#'
#' Liang K-Y, Zeger SL. Longitudinal data analysis using generalized linear
#' models. Biometrika 1986;73(1):13-22.
#' \url{http://doi.org/10.2307/2336267}
#'
#' Halekoh U, Hojsgaard S, Yan, J. The R package geepack for generalized
#' estimating equations. Journal of Statistical Software 2006;15(2):1-11.
#' \url{http://doi.org/10.18637/jss.v015.i02}
#'
#' Kim J-H. Estimating classification error rate: Repeated cross-validation, repeated
#' hold-out and bootstrap. Computational Statistics and Data Analysis
#' 2009:53(11):3735-3745. \url{http://doi.org/10.1016/j.csda.2009.04.009}
#'
#' Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez J-C,
#' Muller M. \code{pROC}: An open-source package for \code{R} and S+ to
#' analyze and compare ROC curves. BMC Bioinformatics 2011;12(77):1-8.
#' \url{http://doi.org/10.1186/1471-2105-12-77}
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' data(unicorns)
#' ## Add a contrived cluster identifier (25 clusters) for demonstration only:
#' uniclus <- unicorns %>%
#'    mutate(cluster = sample(1:25, size = dim(unicorns)[1], replace = TRUE))
#' ## Use gee_screenr:
#' uniobj3 <- gee_screenr(testresult ~ Q1 + Q2 + Q3 + Q5 + Q6 + Q7, id = cluster,
#'                        data = uniclus, link = "logit", Nfolds = 10)
#' class(uniobj3)
#' methods(class = class(uniobj3)[1])
#' methods(class = class(uniobj3)[2])
#' summary(uniobj3)
#' }
#' @import pROC
#' @import geepack
#' @importFrom stats model.frame complete.cases model.response binomial
#' @export
gee_screenr <- function(formula,
                        id = NULL,
                        data = NULL,
                        link = c("logit", "cloglog", "probit"),
                        corstr = c("independence", "exchangeable", "unstructured"),
                        Nfolds = 10,
                        partial_auc = c(0.8, 1.0),
                        partial_auc_focus = "sensitivity",
                        partial_auc_correct = TRUE,
                        boot_n = 4000, conf_level = 0.95,
                        seed = Sys.time(),
                        ...){
    if(!inherits(formula, "formula")) stop("Specify an model formula")
    if(!is.data.frame(data)) stop("Provide a data frame")
    link <- match.arg(link)
    corstr <- match.arg(corstr)
    call <- match.call()
    m <- match(c("formula", "data"), names(call), 0L)
    mf <- call[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mfx <- stats::model.frame(formula, data)
    x <- as.matrix(mfx[, -1])
    dat <- eval(mf, parent.frame())
    dat <- cbind(dat, data.frame(assign(as.character(substitute(id)),
                                        data[[as.character(substitute(id))]])))
    names(dat)[dim(dat)[2]] <- as.character(substitute(id))
    dat <- dat[stats::complete.cases(dat), ]
    if(Nfolds > 0.20*dim(dat)[1])
        stop("Nfolds must be < 20% of number of complete observations")
    y <- stats::model.response(mfx, type = "numeric")
    if(!all(y %in% c(0, 1))) stop("Response variable must be binary (0, 1)")
    prev <- mean(y, na.rm = TRUE)
    geecall <- substitute(geepack::geeglm(formula, id = id,
                                          corstr = corstr, data = dat,
                                          family = binomial(link = link), ...))
    geefit <- eval(geecall)
    parmEst <- coef(geefit)[-1]
    if(any(parmEst < 0))
        warning("Some coefficient(s) < 0; associations should be positive.")
    is.roc <- pROC::roc(geefit$y, as.vector(geefit$fitted.values),
                        auc = TRUE, ci = TRUE, of = "auc",
                        conf.level = conf_level,
                        boot.n = boot_n,
                        partial.auc = partial_auc,
                        partial.auc.focus = partial_auc_focus,
                        partial.auc.correct = partial_auc_correct)
    N <- nrow(dat)
    holdouts <- split(sample(1:N), 1:Nfolds)
    cv.results <- data.frame(NULL)
    cv.coef <- data.frame(NULL)
    X_ho <- data.frame(NULL)
    idc <- as.character(substitute(id))
    set.seed(seed)
    for(i in 1:Nfolds){
        xhoj <- data.frame(fold = rep(i, length(holdouts[[i]])),
                           dat[holdouts[[i]],])
        xhoj <- xhoj[ , -which(names(xhoj) %in% c(idc, as.character(formula[2])))]
        X_ho <- rbind(X_ho, xhoj)
        geecalli <- substitute(geepack::geeglm(formula, id =  id, corstr = corstr,
                                               data = dat[-holdouts[[i]], ],
                                               family = binomial(link = link), ...))
        res = eval(geecalli)
        pred.prob <- inverse_link(stats::predict(res,
                                                newdata = dat[holdouts[[i]], ]),
                                  link = link)
        y <- stats::model.response(mfx[holdouts[[i]], ], type = "numeric")
        cv.results <- rbind(cv.results,
                            data.frame(cbind(fold = rep(i, length(pred.prob)),
                                             y = y,
                                             cv.pred.prob = pred.prob)))
        coef_ <- t(as.matrix(coef(res)))
        colnames(coef_)[1] <- "Intercept"
        coef_ <- data.frame(fold = i, coef_)
        cv.coef <- rbind(cv.coef, coef_)
    }
    attr(X_ho, "Description") <- "Hold-out predictors"
    cv.roc <- pROC::roc(cv.results$y, cv.results$cv.pred.prob,
                        auc = TRUE, ci = TRUE, of = "auc",
                        boot.n = boot_n,
                        conf.level = conf_level,
                        partial.auc = partial_auc,
                        partial.auc.focus = partial_auc_focus,
                        partial.auc.correct = partial_auc_correct)
    class(cv.results) <-  c("cv.predictions", "data.frame")
    result <- list(Call = call,
                   formula = formula,
                   Prevalence = prev,
                   ModelFit = geefit,
                   ISroc = is.roc,
                   Nfolds = Nfolds,
                   CVpreds = cv.results,
                   CVroc = cv.roc,
                   CVcoef = cv.coef,
                   X_ho = X_ho)
    class(result) <- c("gee_screenr", "logreg_screenr")
    invisible(result)
}
