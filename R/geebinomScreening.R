################################################################################
##       R PROGRAM: geebinomScreening.R
##
##         PROJECT: R fuctions for HIV screening
##      WRITTEN BY: Steve Gutreuter, CDC/CGH/DGHT/Statistics, Estimation and
##                               Modeling Team
##                  E-mail:  sgutreuter@cdc.gov
##
##      DISCLAIMER: Although this code has been used by the Centers for Disease
##                  Control & Prevention (CDC), no warranty, expressed or
##                  implied, is made by the CDC or the U.S. Government as to the
##                  accuracy and functioning of the code and related program
##                  material nor shall the fact of distribution constitute any
##                  such warranty, and no responsibility is assumed by the CDC
##                  in connection therewith.
##
################################################################################
#' Test-Screening Tool Based on Marginal Estimation from Mixed-Effects Binomial
#' Models
#'
#' Estimation of the marginal (population averaged) binomial model parameters in
#' the presence of cluster-level random effects, and cross-validated performance
#' of test screening.  The results
#' provide information from which to choose a probability threshold above which
#' individual out-of-sample probabilies indicate the need to perform a diagnostic
#' test.  Out-of-sample performance is estimated using \emph{k}-fold cross
#' validation.
#'
#' @param formula an object of class \code{\link[stats]{formula}}  defining the
#' testing outcome and predictor covariates, which is passed to
#' \code{geepack::geeglm()}.
#' @param id a vector (variable) which identifies the sampling clusters.
#' @param corstr a character string specifying the random-effect correlation
#' structure; one of "exchangeable", "independence", "fixed" or "unstructured".
#' @param data  the "training" sample; a data frame containing the testing
#' outcome and predictive covariates to be used for testing screening.  The
#' testing outcome must be binary (0,1) indicating negative and positive test
#' results, respectively, or logical (\verb{TRUE}/\verb{FALSE}).  The covariates
#' are typically binary (0 = no, 1 = yes) responses to questions which may be
#' predictive of the test result, but any numeric or factor covariates can be
#' used.
#' @param link the character-valued name of the link function for binomial
#' regression.  Choices are \code{"logit"} (default), \code{"cloglog"} or
#' \code{"probit"}.
#' @param Nfolds an integer number of folds used for \emph{k}-fold cross
#' validation (default = 40).
#' @param ... additional arguments passsed to \code{geepack::geeglm}.
#'
#' @return An object of class binomscreenr containing the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{ModelFit}}{An object of class \code{\link[geepack]{geeglm}}}
#' \item{\code{Prevalence}}{Prevalence of the test condition in the training
#' sample.}
#' \item{\code{ParmEst}}{A vector containing the binomial regression parameter
#' estimates.}
#' \item{\code{InSamplePerf}}{A data frame containing in-sample
#' (overly-optimistic) sensitivities and specificities.}
#' \item{\code{CrossVal}}{A data frame containing \emph{k}-fold cross-validation
#' results.}
#' \item{\code{CrossValPerf}}{A data frame containing out-of-sample
#' sensitivities and specificities.}
#' }
#'
#' @seealso \code{\link[geepack]{geeglm}}
#'
#' @examples
#' ## Evaluate the performance of screening thresholds based on a mixed-effect
#' ## logisitc model
#'
#' ##data(unicorns)
#' ##help(unicorns)
#' ##unitool <- geebinomScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5,
#' ##id = "clinic", data = unicorns, link = "logit")
#'
#' ##summary(unitool)
#' ##plot(unitool)
#' ##\dontrun{testCounts(unitool)}
#'
#' ## Example implementation of screening based on those results
#' ## Suppose there are new observations (excluding testing) from two previously
#' ## untested unicorns:
#'
#' ##new <- data.frame(ID = c('"Bernie P."', '"Alice D."'), Q1 = c(0, 0),
#' ##Q2 = c(0, 0), Q3 = c(1, 0), Q4 = c(0, 0), Q5 = c(1, 0),
#' ##clinic = factor(c("C-5", "C-15"),
#' ##                         levels = c("C-1", "C-10", "C-11", "C-12", "C-13",
#' ##                                    "C-14", "C-15", "C-16", "C-17,", "C-18",
#' ##                                    "C-19", "C-2", "C-20", "C-3", "C-4",
#' ##                                    "C-5", "C-6", "C-7", "C-8", "C-9")))
#' ##print(new)
#' ## Compute point estimates of their predicted probabilities testing positive:
#' ##inverseLink("logit",
#' ##            as.matrix(cbind(rep(1, nrow(new)), new[, 2:6])) %*%
#' ##                            as.matrix(unitool$ParamEst, ncol = 1))
#' ## or, more simply,
#' ##predict(unitool$ModelFit, newdata = new, type = "response")
#' ## If, for example, \code{p} = 0.025 is chosen as the screening threshold
#' ## (sensitivity and specificity 77\% and 69\%, respectively) then "Bernie P."
#' ## would be offered testing and "Alice D." would not.
#'
#' ## In practice, the computation of the probabilities of positive test results
#' ## among newly observed individuals might be coded outside of R using, say, a
#' ## spreadsheet.
#' @import pROC
#' @importFrom stats update model.frame complete.cases binomial fitted predict
#' model.response
#' @export
geebinomScreening <- function(formula,
                             id = NULL,
                             corstr = "exchangeable",
                             data = NULL,
                             link = "logit",
                             Nfolds = 40L,
                             ...){
    stop("Function geebinomScreening is broken")   ## FUNCTION BROKEN; See below
    if(!inherits(formula, "formula")) stop("Specify an model formula")
    if(!is.data.frame(data)) stop("Provide a data frame")
    if(!link %in% c("logit", "cloglog", "probit")) stop("Invalid link")
    if(!corstr %in% c("exchangeable", "independence", "fixed", "unstructured"))
        stop("Invalid corstr")
    call <- match.call()
    formx <- update(formula, paste("~ . + ", id))
    mf <- model.frame(formx, data = data)
    dat <- eval(mf, parent.frame())
    dat <- dat[complete.cases(dat), ]
    dat <- dat[order(dat[[id]]), ]
    if(Nfolds > 0.20*dim(dat)[1])
        stop("Nfolds must be < 20% of number of complete observations")
    y <- stats::model.response(dat, "numeric")
    if(!all(y %in% c(0, 1))) stop("Response variable must be binary (0, 1)")
    prev <- mean(y, na.rm = TRUE)
    ## FIXME: cls <- model.extract(mf, clinic) fails.  Is mf a model.frame?
    ##cls <- model.extract(mf, clinic)
    lrfit <- geepack::geeglm(formula, id = id, corstr = corstr, data = dat,
                      family = binomial(link = link))
    is.roc <- pROC::roc(lrfit@resp$y, fitted(lrfit))
    N <- nrow(dat)
    holdouts <- split(sample(1:N), 1:Nfolds)
    cv.results <- data.frame(NULL)
    for(i in 1:Nfolds){
        res <- geepack::geeglm(formula, id = dat[[id]], corstr = corstr,
                        data = dat[-holdouts[[i]], ],
                        family = binomial(link = link))
        pred.prob <- inverseLink(link,
                                 predict(res, newdata = dat[holdouts[[i]], ]))
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
                   ParamEst = lrfit@beta,
                   ISroc = is.roc,
                   CVpreds = cv.results,
                   CVroc = cv.roc)
    class(result) <- "binomscreenr"
    invisible(result)
}


################################   END of FILE   ################################
