#################################################################################
##       R PROGRAM: binomialScreening.R
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

#' Development of a Testing Screening Tool Based on Binomial Regression
#'
#' \code{binomialScreening} Development of a screeing tool based on binomial
#' regression.
#'
#' Development of a screening tool based on predicted probabilities of positive
#' test results based on logistic regression.  The results provide information
#' from which to choose a probability threshold above which individual
#' out-of-sample predictive probabilies indicate the need to test that
#' individual.  Out-of-sample performance is estimated using \emph{k}-fold
#' cross validation.
#'
#' S3 \code{plot}, \code{print} and \code{summary} methods are available.
#'
#' @param formula A formula which is passed to \code{stats::glm()}.
#' @param data  A data frame containing the testing outcome and predictive
#' covariates to be used for testing screening.  The testing outcome must
#' be binary (0,1) indicating negative and positive test results, respectively.
#' The covariates are typically binary (0 = no, 1 = yes) responses to
#' questions which may be predictive of the test result, but any covariates
#' compatible with \code{glm} can be used.
#' @parm link Link function for binomial regression.  Choices are "logit"
#' (default), "cloglog" or "probit".
#' @param Nfolds The number of folds used for \emph{k}-fold cross validation.
#' Default = 10.
#' @param p.threshold A vector of reference probabilities for estimation of
#' Receiver Operating Characteristics.
#'
#' @return A list of class LRscreeing containing the elements:
#' \describe{
#' \item{Call} The function call.
#' \item{ModelFit} An object produced by \code{stats::glm()}.
#' \item{ParmEst} The binomial regression  parameter estimates.
#' \item{InSamplePerf} A data frame containing in-sample sensitivities and
#' specificities.
#' \item{CrossVal} A data frame containing cross-validation results.
#' \item{CrossValPerf} A data frame containing out-of-sample
#' sensitivities and specificities.
#' }
#'
#' @author Steve Gutreuter, \email{sgutreuter@@gmail.com}
#'
#' @seealso \code{\link{stats::glm}}
#'
#' @examples
#' data(unicorns)
#' unitool <- binomialScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5,
#'                              data = unicorns, Nfolds = 20,
#'                              p.threshold = c(seq(0.01, 0.10, by = 0.015),
#'                                              seq(.15, 0.95, by = 0.05)))
#' summary(unitool)
#' plot(wtf)
#' @export
binomialScreening <- function(formula, data = NULL, link = "logit", Nfolds = 10,
                              p.threshold = c(seq(0.01, 0.09, 0.01),
                                              seq(0.1, 0.6, 0.05),
                                              seq(0.65, 0.95, 0.10 ))){
    if(!plyr::is.formula(formula)) stop("Specify an model formula")
    if(!is.data.frame(data)) stop("Provide a data frame")
    if(!link %in% c("logit", "cloglog", "probit")) stop("Invalid link")
    call <- match.call()
    m <- match(c("formula", "data"), names(call), 0L)
    mf <- call[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    dat <- eval(mf, parent.frame())
    dat <- dat[complete.cases(dat), ]
    if(Nfolds > 0.20*dim(dat)[1])
        stop("Nfolds must be < 20% of number of complete observations")
    y <- model.response(dat, "numeric")
    if(!all(y %in% c(0, 1))) stop("Response variable must be binary (0, 1)")
    lrfit <- stats::glm(formula, data = dat, family = binomial(link = link))
    insamp <- data.frame(NULL)
    pp <- inverseLink(link, lrfit$linear.predictors)
    for(j in seq_along(p.threshold)){
        I.test <-as.numeric(pp >= p.threshold[j])
        cmat <- table(factor(I.test, levels = c("0", "1")),
                      factor(y, levels = c("0", "1")))
        z <- sens_spec(cmat)
        insamp <- rbind(insamp, data.frame(p.threshold = p.threshold[j],
                                               sensitivity = z[1],
                                               specificity = z[2]))
    }
    cv.results <- data.frame(NULL)
    N <- nrow(dat)
    holdouts <- split(sample(1:N), 1:Nfolds)
    for(i in 1:Nfolds){
        res <- stats::glm(formula, data = dat[-holdouts[[i]], ],
                          family = binomial(link = link))
        pred.prob <- inverseLink(link, predict(res, newdata = dat[holdouts[[i]], ]))
        y <- model.response(dat[holdouts[[i]], ])
        cv.results <- rbind(cv.results,
                            data.frame(cbind(fold = rep(i, length(pred.prob)),
                                             y = y,
                                             data = dat[holdouts[[i]],],
                                             cv.pred.prob = pred.prob)))
    }
    SensSpec.results <- data.frame(NULL)
    for(j in seq_along(p.threshold)){
        I.test <-as.numeric(cv.results$cv.pred.prob >= p.threshold[j])
        cmat <- table(factor(I.test, levels = c("0", "1")),
                      factor(cv.results$y, levels = c("0", "1")))
        z <- sens_spec(cmat)
        SensSpec.results <- rbind(SensSpec.results,
                                  data.frame(p.threshold = p.threshold[j],
                                             sensitivity = z[1],
                                             specificity = z[2]))
    }
    class(insamp) <-  c("ROC", "data.frame")
    class(SensSpec.results) <-  c("ROC", "data.frame")
    result <- list(Call = call,
                   ModelFit = lrfit,
                   ParmEst = lrfit$coeff,
                   InSamplePerf = insamp,
                   CrossVal = cv.results,
                   CrossValPerf = SensSpec.results)
    class(result) <- c("binomscreenr", "list")
    result
}



#' An S3 Summary Method
#'
#' \code{summary.binomscreenr} A \code{summary} method for \code{binomscreenr}
#' objects.
#'
#' @param obj An object of class \code{binomscreenr} produced by function
#' \code{binomialScreening}
#'
#' @author Steve Gutreuter, \email{sgutreuter@@gmail.com}
#'
#' @examples
#' @export
summary.binomscreenr <- function(obj){
    if(!("binomscreenr" %in% class(obj))) stop("obj not binomscreenr class")
    cat("Call:\n")
    print(obj$Call)
    cat("\n\nLogistic regression model summary:")
    print(summary(obj$ModelFit))
    cat("\nOut-of-sample sensitivity and specificity\nscreening at p.hat >= p.threshold:\n\n")
    print(obj$CrossValPerf)
}
#' An S3 Print Method
#'
#' \code{print.binomscreenr} A print method for \code{binomscreenr} objects.
#'
#' @param obj An object of class \code{binomscreenr} produced by function
#' \code{binomialScreening}
#'
#' @author Steve Gutreuter, \email{sgutreuter@@gmail.com}
#'
print.binomscreenr <- function(obj){
    if(!("binomscreenr" %in% class(obj))) stop("obj not binomscreenr class")
    cat("Out-of-sample sensitivity and specificity\nscreening at p.hat >= p.threshold:\n\n")
    print(obj$CrossValPerf)
}


#' An S3 Plot Method
#'
#' \code(plot.binomscreenr) A \code{plot} method for \code{binomscreenr}
#' objects.
#'
#' A wrapper function for \code{lattice::xyplot} customized to produce a plot
#' of Receiver Operating Characteristic curves.
#'
#' @param obj An object of class \code{binomscreenr} produced by function
#' \code{binomialScreening}
#'
#' @return A \code{lattice} graphical object
#'
#' @author Steve Gutreuter, \email{sgutreuter@@gmail.com}
#'
#' @export
plot.binomscreenr <- function(obj, main = "Receiver Operating Characteristics",
                             ...){
    nrows <- dim(obj$CrossValPerf)[1]
    d.lty <- c(2, 1)
    dfrm <- data.frame(p.threshold =
                           rep(obj$CrossValPerf$p.threshold, 2),
                       grp = c(rep("In-sample", nrows),
                                   rep("Out-of-sample", nrows)),
                       sensitivity  = c(obj$InSamplePerf$sensitivity,
                                    obj$CrossValPerf$sensitivity),
                       FPP = 1 - c(obj$InSamplePerf$specificity,
                                        obj$CrossValPerf$specificity))
    d.key <- list(corner = c(1, 0), x = 0.98, y = 0.04,
                  text = list(c("In-sample (overly-optimistic)",
                                "Out-of-sample")),
                  lines = list(type = c("l", "l"), col = rep("black", 2),
                               lty = d.lty))
    res <- lattice::xyplot(sensitivity ~ FPP,
                           group = grp,
                           data = dfrm,
                           main = main,
                           type = rep("S", 2),
                           col = rep("black", 2),
                           lty = d.lty,
                           ylab = "Sensitivity (%)",
                           xlab = "1 - Specificity (%)",
                           key = d.key)
    res
}

################################   END of FILE   ################################
