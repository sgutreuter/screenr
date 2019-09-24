#################################################################################
##       R PROGRAM: mebinomScreening.R
##
##         PROJECT: R fuctions for HIV screening tool development
##
##      WRITTEN BY: Steve Gutreuter, CDC/CGH/DGHT/Statistics, Estimation and
##                               Modeling Team
##                  E-mail:  sgutreuter@cdc.gov
##
##      DISCLAIMER: Although this code has been used by the Centers for Disease
##                  Control & Prevention (CDC), no warranty, expressed or
##                  implied, is made by the CDC or the U.S. Government as to the
##                  accuracy and functioning of the code and related program
##                  material nor shall the fact of distribution constitute any
##                  such warranty, and no responsibility is assumed by the CDC in
##                  connection therewith.
##
#################################################################################
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
#' @param formula an object of class \code{\link{formula}}  defining the testing
#' outcome and predictor covariates, which is passed to \code{lme4::glmer()}.
#' @param id a vector (variable) which identifies the sampling clusters.
#' @param data  the "training" sample; a data frame containing the testing
#' outcome and predictive covariates to be used for testing screening.  The
#' testing outcome must be binary (0,1) indicating negative and positive test
#' results, respectively, or logical (TRUE/FALSE).  The covariates are typically
#' binary (0 = no, 1 = yes) responses to questions which may be predictive of
#' the test result, but any numeric or factor covariates can be used.
#' @param link the character-valued name of the link function for binomial
#' regression.  Choices are "\code{logit}" (default), "\code{cloglog}" or
#' "\code{probit}".
#' @param Nfolds an integer number of folds used for \emph{k}-fold cross
#' validation (default = 20).
#' @param p a numeric vector of reference probabilities for estimation of
#' Receiver Operating Characteristics.
#' @param ... additional arguments passsed to \code{lme4::glmer}.
#'
#' @return An object of class binomscreenr containing the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{ModelFit}}{An object of class \code{\link{glmerMod}}}
#' \item{\code{Prevalence}}{Prevalence of the test condition in the training sample.}
#' \item{\code{ParmEst}}{A vector containing the binomial regression parameter estimates.}
#' \item{\code{InSamplePerf}}{A data frame containing in-sample (overly-optimistic)
#' sensitivities and specificities.}
#' \item{\code{CrossVal}}{A data frame containing \emph{k}-fold cross-validation results.}
#' \item{\code{CrossValPerf}}{A data frame containing out-of-sample  sensitivities and
#' specificities.}
#' }
#'
#' @seealso \code{\link{lme4::glmer}}
#'
#' @references
#'
#'
#' @examples
#' ## Evaluate the performance of screening thresholds based on a logisitc model
#'
#' data(unicorns)
#' help(unicorns)
#' unitool <- mebinomScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5,
#'                              data = unicorns, link = "logit",
#'                              p = c(seq(0.01, 0.10, by = 0.015),
#'                                    seq(.15, 0.95, by = 0.05)))
#' summary(unitool)
#' plotROC(unitool)
#' testCounts(SensSpec = unitool)
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
#' ## Compute point estimates of their predicted probabilities testing positive:
#' inverseLink("logit",
#'             as.matrix(cbind(c(nrow(new), nrow(new)), new[, 2:6])) %*%
#'                             as.matrix(unitool$ParmEst, ncol = 1))
#' ## or, more simply,
#' predict(unitool$ModelFit, newdata = new, type = "response")
#' ## If, for example, \code{p} = 0.025 is chosen as the screening threshold
#' ## (sensitivity and specificity 77\% and 69\%, respectively) then "Bernie P."
#' ## would be offered testing and "Alice D." would not.
#'
#' ## In practice, the computation of the probabilities of positive test results
#' ## among newly observed individuals might be coded outside of R using, say, a
#' ## spreadsheet.
#'
#' @export
mebinomScreening <- function(formula,
                         id = NULL,
                         data = NULL,
                         link = "logit",
                         Nfolds = 20L,
                         p = c(seq(0.01, 0.09, 0.01),
                               seq(0.1, 0.6, 0.05),
                               seq(0.65, 0.95, 0.10 )),
                         ...){
    ##stop("This function is temporarily out of order.")   ###############################
    if(!plyr::is.formula(formula)) stop("Specify an model formula")
    if(!is.data.frame(data)) stop("Provide a data frame")
    if(!link %in% c("logit", "cloglog", "probit")) stop("Invalid link")
    call <- match.call()
    meform <- update(formula, paste("~ . + (1|", id, ")"))
    formx <- update(formula, paste("~ . +", id))
    mf <- model.frame(formx, data)
    ## m <- match(c("formula", "data"), names(call), 0L)
    ## mf <- call[c(1L, m)]
    ## mf[[1L]] <- quote(stats::model.frame)
    dat <- eval(mf, parent.frame())
    cc <- complete.cases(dat)
    dat <- dat[cc, ]
    if(Nfolds > 0.20*dim(dat)[1])
        stop("Nfolds must be < 20% of number of complete observations")
    y <- stats::model.response(mf, "numeric")
    if(!all(y %in% c(0, 1))) stop("Response variable must be binary (0, 1)")
    prev <- mean(y, na.rm = TRUE)
    modfit <- lme4::glmer(meform, data = dat, family = binomial(link = link))
    insamp <- data.frame(NULL)
    browser()                          #########################################################
    pp <- fitted(modfit)               ############   Check this   #############################
    pp2 <- predict(modfit, type = "response")   ######### remove ###############################
    View(cbind(pp, pp2))  ####### remove ##############################
    for(j in seq_along(p)){
        I.test <-as.numeric(pp >= p[j])
        cmat <- table(factor(I.test, levels = c("0", "1")),
                      factor(y, levels = c("0", "1")))
        z <- sens_spec(cmat)
        insamp <- rbind(insamp, data.frame(p = p[j],
                                           sensitivity = z[1],
                                           specificity = z[2]))
    }
    cv.results <- data.frame(NULL)
    N <- nrow(dat)
    holdouts <- split(sample(1:N), 1:Nfolds)
    for(i in 1:Nfolds){
        dati <- dat[-holdouts[[i]], ]
        res <- lme4::glmer(meform, data = dati, family = binomial(link = link))
        pred.prob <- predict(res, newdata = dat[holdouts[[i]], ], re.form = NULL,
                             type = "response")
    y <- model.response(dat[holdouts[[i]], ])
        cv.results <- rbind(cv.results,
                            data.frame(cbind(fold = rep(i, length(pred.prob)),
                                             y = y,
                                             data = dat[holdouts[[i]],],
                                             cv.pred.prob = pred.prob)))
    }
    SensSpec.results <- data.frame(NULL)
    for(j in seq_along(p)){
        I.test <-as.numeric(cv.results$cv.pred.prob >= p[j])
        cmat <- table(factor(I.test, levels = c("0", "1")),
                      factor(cv.results$y, levels = c("0", "1")))
        z <- sens_spec(cmat)
        SensSpec.results <- rbind(SensSpec.results,
                                  data.frame(p = p[j],
                                             sensitivity = z[1],
                                             specificity = z[2]))
    }
    class(insamp) <-  c("ROC", "data.frame")
    class(SensSpec.results) <-  c("ROC", "data.frame")
    result <- list(Call = call,
                   ModelFit = modfit,
                   Prevalence = prev,
                   ParmEst = modfit@beta,
                   InSamplePerf = insamp,
                   CrossVal = cv.results,
                   CrossValPerf = SensSpec.results)
    class(result) <- "binomscreenr"
    result
}



#' An S3 Summary Method for \code{binomscreenr} Objects.
#'
#' @param object an object of class \code{binomscreenr} produced by function
#' \code{binomialScreening}
#' @param ... further arguments passed to or from other methods.
#'
#' @export
summary.binomscreenr <- function(object, ...){
    if(!("binomscreenr" %in% class(object))) stop("object not binomscreenr class")
    cat("Call:\n")
    print(object$Call)
    cat("\n\nLogistic regression model summary:")
    print(summary(object$ModelFit))
    cat("\nOut-of-sample sensitivity and specificity\nscreening at p.hat >= p:\n\n")
    print(object$CrossValPerf)
}
#' An S3 Print Method for \code{binomscreenr} Objects.
#'
#' @param x an object of class \code{binomscreenr} produced by function.
#' @param ... further arguments passed to or from other methods.
#' @param quote logical, indicating whether or not strings should be printed
#' with surrounding quotes.
#' @export
print.binomscreenr <- function(x, quote = FALSE, ...){
    if(!("binomscreenr" %in% class(x))) stop("x not binomscreenr class")
    cat("Out-of-sample sensitivity and specificity\nscreening at p.hat >= p:\n\n")
    print(x$CrossValPerf)
}


################################   END of FILE   ################################
