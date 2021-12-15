#################################################################################
##       R PROGRAM: logreg_screenr.R
##
##         PROJECT: R function for HIV screening based on ordinary logistic
##                  regression.
##
##      WRITTEN BY: Steve Gutreuter
##                  E-mail:  sgutreuter@gmail.com
##
#################################################################################

## Function logreg_screenr
##
#' Fitting Screening Tools Using Ordinary Logistic Regression
#'
#' @description
#' \code{logreg_screenr} is a convenience function which integrates ordinary
#' logistic regression \emph{k}-fold  cross-validation and estimation of the
#' receiver-operating characteristic.
#'
#' @param formula an object of class \code{stats::formula} defining the
#' testing outcome and predictor covariates, which is passed to \code{stats::glm()}.
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
#' @param seed random-number generator seed for cross-validation data splitting.
#'
#' @param ... additional arguments passsed to or from other \code{stats::glm}
#' or \code{pROC::roc}.
#'
#' @return An object of class logreg_screenr containing the elements:
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
#' For a gentle but python-centric introduction to \emph{k}-fold cross-validation,
#' see \url{https://machinelearningmastery.com/k-fold-cross-validation/}.
#'
#' @note
#' \code{logreg_screenr} is intended mainly for comparison with \code{lasso_screenr}.
#' Careful manual model selection is required with \code{logreg_screenr}.
#' \code{lasso_screenr} is easier and should generally produce better results.
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
#' uniobj2 <- logreg_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
#'                            data = unicorns, link = "logit", Nfolds = 10)
#' summary(uniobj2)
#'
#' @aliases{binomialScreenr}
#' @import pROC
#' @importFrom stats model.frame complete.cases model.response glm binomial
#' @export
logreg_screenr <- function(formula,
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
    mfx <- stats::model.frame(formula, data)
    x <- as.matrix(mfx[, -1])
    dat <- eval(mf, parent.frame())
    dat <- dat[stats::complete.cases(dat), ]
    if(Nfolds > 0.20*dim(dat)[1])
        stop("Nfolds must be < 20% of number of complete observations")
    y <- stats::model.response(dat, "numeric")
    if(!all(y %in% c(0, 1))) stop("Response variable must be binary (0, 1)")
    prev <- mean(y, na.rm = TRUE)
    lrfit <- stats::glm(formula, data = dat, family = binomial(link = link), ...)
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
                          family = binomial(link = link), ...)
        pred.prob <- inverse_link(stats::predict(res,
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
    class(result) <- "logreg_screenr"
    invisible(result)
}
