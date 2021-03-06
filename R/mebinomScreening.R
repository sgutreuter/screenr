#################################################################################
##       R PROGRAM: mebinomScreening.R
##
##         PROJECT: R fuctions for HIV screening
##      WRITTEN BY: Steve Gutreuter
##                  E-mail:  sgutreuter@gmail.com
##
#################################################################################
#' Test-Screening Tool Based on Marginal Estimation from Mixed-Effects Binomial
#' Models
#'
#' Estimation of the mixed-effect binomial model parameters in
#' the presence of cluster-level random effects, and cross-validated performance
#' of test screening.  The approximate "population-averaged" parameter estimates
#' are obtained by ignoring the random effects after model fitting.
#'
#' The results
#' provide information from which to choose a probability threshold above which
#' individual out-of-sample probabilies indicate the need to perform a diagnostic
#' test.  Out-of-sample performance is estimated using \emph{k}-fold cross
#' validation.
#'
#' The receiver operating characteristics are computed using the \code{pROC}
#' package. See References and package documentation for additional details.
#'
#' @param formula an object of class \code{\link[stats]{formula}}  defining the
#' testing outcome and predictor covariates, which is passed to
#' \code{lme4::glmer()}.
#' @param id a character-valued name of the variable which identifies the
#' sampling clusters.
#' @param data  the "training" sample; a data frame containing the testing
#' outcome and predictive covariates to be used for testing screening.  The
#' testing outcome must be binary (0,1) indicating negative and positive test
#' results, respectively, or logical (\verb{TRUE}/\verb{FALSE}).  The covariates
#' are typically binary (0 = no, 1 = yes) responses to questions which may be
#' predictive of the test result, but any numeric or factor covariates can be
#' used.
#' @param link the character-valued name of the link function for binomial
#' regression.  Choices are \verb{"logit"} (default), \verb{"cloglog"} or
#' \code{"probit"}.
#' @param Nfolds number of folds used for \emph{k}-fold cross
#' validation (default = 40).
#' @param ... additional arguments passsed to \code{lme4::glmer}.
#'
#' @return An object of class \code{binomscreenr} containing the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{ModelFit}}{An object of class \code{\link[lme4]{merMod-class}}}
#' \item{\code{Prevalence}}{Prevalence of the test condition in the training sample.}
#' \item{\code{ParmEst}}{A vector containing the binomial regression parameter estimates.}
#' \item{\code{ISroc}}{An object of class \code{\link[pROC]{roc}} containing
#' the "in-sample" (overly-optimistic) receiver operating characteristics,
#' and additional functions for use with this object are
#' available in the \code{pROC} package.}
#' \item{\code{CVpreds}}{An object of class \code{cv.predictions} containing
#' the data and cross-validated predicted condition \code{y}.}
#' \item{\code{CVroc}}{An object of class \code{\link[pROC]{roc}} containing
#' the \emph{k}-fold cross-validated "out-of-sample" receiver operating
#' characteristics, and additional functions for use with this object are
#' available in the \code{pROC} package.}
#' }
#'
#' @seealso \code{\link[lme4]{glmer}}
#'
#' @references
#' Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez J-C,
#' Müller M. \code{pROC}: An open-source package for \code{R} and S+ to
#' analyze and compare ROC curves. BMC Bioinformatics. 2011;12(77):1-8.
#' \url{http://doi.org/10.1186/1471-2105-12-77}
#'
#' @examples
#' ## Evaluate the performance of screening thresholds based on a mixed-effect
#' ## logisitc model
#'
#' data(unicorns)
#' help(unicorns)
#' unitool <- mebinomScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5, id = "clinic",
#'                              data = unicorns, link = "logit", Nfolds = 5L)
#' summary(unitool)
#' plot(unitool)
#' \dontrun{testCounts(unitool)}
#'
#' ## Example implementation of screening based on those results
#' ## Suppose there are new observations (excluding testing) from two previously
#' ## untested unicorns:
#'
#' new <- data.frame(ID = c('"Bernie P."', '"Alice D."'), Q1 = c(0, 0), Q2 = c(0, 0),
#'                    Q3 = c(1, 0), Q4 = c(0, 0), Q5 = c(1, 0),
#'                    clinic = factor(c("C-5", "C-15"),
#'                                    levels = c("C-1", "C-10", "C-11", "C-12", "C-13",
#'                                               "C-14", "C-15", "C-16", "C-17,", "C-18",
#'                                               "C-19", "C-2", "C-20", "C-3", "C-4",
#'                                               "C-5", "C-6", "C-7", "C-8", "C-9")))
#' print(new)
#' ## Compute the predicted probabilities that the new unicorns will test
#' ## positive, ignoring the random effects:
#' predict(unitool$ModelFit, newdata = new, type = "response", re.form = ~0)
#' ## If, for example, \code{p} = 0.025 is chosen as the screening threshold
#' ## (sensitivity and specificity 77\% and 69\%, respectively) then "Bernie P."
#' ## would be offered testing and "Alice D." would not.
#'
#' ## In practice, the computation of the probabilities of positive test results
#' ## among newly observed individuals might be coded outside of R using, say, a
#' ## spreadsheet.
#' @import lme4 pROC
#' @importFrom stats update model.frame complete.cases binomial fitted predict model.response
#' @export
mebinomScreening <- function(formula,
                             id = NULL,
                             data = NULL,
                             link = "logit",
                             Nfolds = 40L,
                              ...){
    if(!inherits(formula, "formula")) stop("Specify an model formula")
    if(!is.data.frame(data)) stop("Provide a data frame")
    if(!link %in% c("logit", "cloglog", "probit")) stop("Invalid link")
    call <- match.call()
    meform <- update(formula, paste("~ . + (1|", id, ")"))
    formx <- update(formula, paste("~ . +", id))
    mf <- stats::model.frame(formx, data)
    dat <- eval(mf, parent.frame())
    dat <- dat[complete.cases(dat), ]
    if(Nfolds > 0.20*dim(dat)[1])
        stop("Nfolds must be < 20% of number of complete observations")
    y <- stats::model.response(dat, "numeric")
    if(!all(y %in% c(0, 1))) stop("Response variable must be binary (0, 1)")
    prev <- mean(y, na.rm = TRUE)
    lrfit <- lme4::glmer(meform, data = dat, family = binomial(link = link))
    is.roc <- pROC::roc(lrfit@resp$y, fitted(lrfit))
    N <- nrow(dat)
    holdouts <- split(sample(1:N), 1:Nfolds)
    cv.results <- data.frame(NULL)
    for(i in 1:Nfolds){
        res <- lme4::glmer(meform, data = dat[-holdouts[[i]], ],
                           family = binomial(link = link))
        pred.prob <- predict(res,
                             newdata = dat[holdouts[[i]], ],
                             type = "response",
                             re.form = ~0)
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
                   ParamEst = lrfit@beta,
                   ISroc = is.roc,
                   CVpreds = cv.results,
                   CVroc = cv.roc)
    class(result) <- "binomscreenr"
    invisible(result)
}


################################   END of FILE   ################################
