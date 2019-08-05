#################################################################################
##       R PROGRAM: testScreening.R
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


#' Sample Size for Joint Testing of Sensitivity and Specificity
#'
#' \code{nSensSpec} Returns required number of positives and negatives
#'
#' Sample size required for testing against the joint null hypothesis
#' H0: Sens = SnsCrit and Spec = SpcCrit
#'
#' @param Sens A vector of snticipated sensitivities (0 < Sens < 1).
#' @param Spec A vector of anticipated specificities (0 < Sens < 1).
#' @param SnsCrit Critical value for sensitivity.
#' @param SnsCrit Critical value for specificity.
#' @param alpha Probability of Type I error.
#' @param power Desired power.
#' @return A data frame containing the required number of positive and negative
#' test results.
#' @references Sullivan, M.P., The Statistical Evaluation of Medical Tests
#' for Classification and Prediction.  Oxford Statistical Science Series 28.
#' Oxford University Press, Oxford UK.
#' @examples
#' sens <- c(rep(0.8, 4), rep(0.85, 4), rep(0.90, 4))
#' spec <- rep(c(0.65, 0.70, 0.75, 0.8 ), 3)
#' nSensSpec(sens, spec, SnsCrit = 0.70, SpcCrit = 0.60)
#' @export
nSensSpec <- function(Sens, Spec, SnsCrit = 0.9, SpcCrit = 0.9,
                      alpha = 0.05, power =0.8){
    if(any(Sens <0 | Sens > 1)) stop("Sensitivity not in [0,1]")
    if(any(Spec <0 | Spec > 1)) stop("Specificity not in [0,1]")
    alpha.star <- 1 - sqrt(1 - alpha)
    beta.star <- 1 - sqrt(power)
    Z <- qnorm(c((1 - alpha.star), (1 - beta.star)))
    n.pos <- (Z[1]*sqrt(SnsCrit*(1 - SnsCrit)) + Z[2]*sqrt(Sens*(1 - Sens)))^2 /
        (Sens - SnsCrit)^2
    n.neg <- (Z[1]*sqrt(SpcCrit*(1 - SpcCrit)) + Z[2]*sqrt(Spec*(1 - Spec)))^2 /
        (Spec - SpcCrit)^2
    data.frame(n.pos = ceiling(n.pos), n.neg = ceiling(n.neg))
}

#' Sensitivity and Specificity from a 2 x 2 Table
#'
#' \code{SensSpec} Returns sensitivity and specificity of a test.
#' @parm x A 2 x 2 table, with negatives appearing first in rows and columns.
#' @return A list containing components sensitivity and specificity.
#' @examples
#' TrueResults <- ordered(c(0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0))
#' TestResults <- ordered(c(0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0))
#' sens_spec(table(TestResults, TrueResults))
#' @export
sens_spec <- function(x){
    if(!class(x) == "table") stop('arg class is not "table"')
    if(!(dim(x)[1] == 2L & dim(x)[2] == 2L)) stop("arg not a 2x2 table")
    spec = x[1, 1] / sum(x[ , 1])
    sens = x[2, 2] / sum(x[ , 2])
    result = list(sensitivity = sens, specificity = spec)
    result
}

#' A Simple Un-optimized HIV Testing Screener
#'
#' \code{Bandason} Development of a simple screening tool.
#'
#' Development of a simple testing screening tool based on the sum of ordinal
#' numeric responses to questions which are predictive of the testing results.
#' \code{logisticScreening} will generally outperform this approach, and
#' \code{Bandason} is provided mainly for comparison.
#'
#' S3 \code{plot}, \code{print} and \code{summary} methods are available.
#'
#' @param formula A formula for the screener.
#' @param data A data frame containing the testing outcome and predictive
#' covariates to be used for testing screening.  The testing outcome must
#' be binary (0,1) indicating negative and positive test results, respectively.
#' The covariates are typically binary (0 = no, 1 = yes) responses to
#' questions which may be predictive of the test result.
#'
#' @return A list having S3 class "Bandason" containing the call, the
#' Receiver Operating Characteristics and the scores.
#'
#' @author Steve Gutreuter, \email{sgutreuter@@gmail.com}
#'
#' @references Bandason et al 2016. Validation of a screening tool to
#' identify older children living with HIV in primary care facilities in high
#' HIV prevalence settings. AIDS 30(5):779-785
#' \url{http://dx.doi.org/10.1097/QAD.0000000000000959
#'
#' @examples
#' TO DO...
#'
#' @export
Bandason <- function(formula, data){
    mf <- match.call(expand.dots = FALSE)
    call <- mf
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    dat <- data[complete.cases(data[, names(mf), ]), ]
    y <- model.response(mf, "numeric")
    preds <- mf[, -1]
    npreds <- dim(preds)[2]
    score <- apply(preds, 1, sum)
    critvals <- 1L:npreds
    tst <- outer(score, critvals, `>=`)
    roc <- data.frame(score = NULL, sensitivity = NULL, specificity = NULL)
    for(i in seq_along(critvals)){
        tabl <- table(factor(tst[, i] * 1 , levels = c("0", "1")),
                      factor(y, levels = c("0", "1")))
        ss <- sens_spec(tabl)
        roc <- rbind(roc, data.frame(score = i,
                                     sensitivity = ss$sensitivity,
                                     specificity = ss$specificity))
        class(roc) <- c("ROC", "data.frame")
    }
    scores <- cbind(dat, Bandason.score = score)
    result <- list(Call = call, roc = roc, scores = scores)
    class(result) <- "Bandason"
    invisible(result)
}


#' An S3 Print Method
#'
#' \code(print.Bandason) An S3 printing method for objects created by
#' function \code(Bandason()).
#'
#' @param x A Bandason class object.
#'
#' @examples
#' TO DO
print.Bandason <- function(x){
    if(!class(x) == "Bandason") stop("x not a Bandason object")
    x$roc
}

#' An S3 Plot Method
#'
#' \code(plot.Bandason) An S3 plotting method for objects created by
#' function \code(Bandason()).
#'
#' @param x A Bandason class object.
#'
#' @examples
#' TO DO
plot.Bandason <- function(x){
    if(!("Bandason" %in% class(x))) stop("x must be a Bandason object")
    plot.ROC(x$roc)
}

#' An S3 Summary Method
#'
#' \code{summary.Bandason} An S3 summary method for objects created by
#' function \code{Bandason()}.
#'
#' @param x A Bandason class object.
#'
#' @examples
#' TO DO
summary.Bandason <- function(x){
    if(!class(x) == "Bandason") stop("x not a Bandason object")
    cat("A list object with three components:\n")
    cat("\nCall:\n")
    print(x$Call)
    cat("\nroc (In-sample Receiver Operating Characteristic):\n")
    print(x$roc)
    cat("\nscores (data frame including Bandason scores)\n")
}

#' #' An S3 Plot Method
#'
#' \code{plot.ROC} An S3 plotting method for Receiver Operating Characteristic
#' objects created by function plot.Bandason().
#'
#' @param x A Bandason class object.
#'
#' @examples
#' TO DO
plot.ROC <- function(x){
    if(!("ROC" %in% class(x))) stop("x must be an ROC object")
    x$FPP <- 1 - x$specificity
    plot(x$FPP, x$sensitivity, type = "S", ylim = c(0,1),
         xlim = c(0, 1), ylab = "Sensitivity", xlab = "1 - Specificity",
         main = "Receiver Operating Characteristic")
    abline(a = 0, b = 1, lty = 2)
    text(x$sensitivity ~ x$FPP, labels = x$score, pos = 2)
}


#' Development of a Testing Screening Tool Based on Logistic Regression
#'
#' \code{logisticScreening} Development of a screeing tool based on logistic
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
#' @param formula A formula which is passed to \code{glm()}.
#' @param data  A data frame containing the testing outcome and predictive
#' covariates to be used for testing screening.  The testing outcome must
#' be binary (0,1) indicating negative and positive test results, respectively.
#' The covariates are typically binary (0 = no, 1 = yes) responses to
#' questions which may be predictive of the test result, but any covariates
#' compatible with \code{glm} can be used.
#' @param Nfolds The number of folds used for \emph{k}-fold cross validation.
#' Default = 10.
#' @param p.threshold A vector of reference probabilities for estimation of
#' Receiver Operating Characteristics.
#'
#' @return A list of class LRscreeing containing the elements:
#' \describe{
#' \item{Call} The function call.
#' \item{ModelFit} An object produced by \code{glm()}.  Any methods for
#' \code{glm} may be used for this element.
#' \item{ParmEst} The logistic parameter estimates.
#' \item{InSamplePerf} A data frame containing in-sample sensitivities and
#' specificities.
#' \item{CrossVal} A data frame containing cross-validation results.
#' \item{CrossValPerf} A data frame containing out-of-sample
#' sensitivities and specificities.
#' }
#'
#' @author Steve Gutreuter, \email{sgutreuter@@gmail.com}
#'
#' @seealso \code{\link{glm}}
#'
#' @examples
#' TO DO
#' @export
logisticScreening <- function(formula, data = NULL, Nfolds = 10,
                              p.threshold = c(seq(0.01, 0.09, 0.01),
                                              seq(0.1, 0.6, 0.05),
                                              seq(0.65, 0.95, 0.10 ))){
    if(!plyr::is.formula(formula)) stop("Specify an model formula")
    if(!is.data.frame(data)) stop("Provide a data frame")
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
    lrfit <- glm(formula, data = dat, family = "binomial")
    insamp <- data.frame(NULL)
    pp <- locfit::expit(lrfit$linear.predictors)
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
        res <- glm(formula, data = dat[-holdouts[[i]], ], family = "binomial")
        pred.prob <- locfit::expit(predict(res, newdata = dat[holdouts[[i]], ]))
        ## y goes bad here; need y from holdouts
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
    result <- list(Call = call,
                   ModelFit = lrfit,
                   ParmEst = lrfit$coeff,
                   InSamplePerf = insamp,
                   CrossVal = cv.results,
                   CrossValPerf = SensSpec.results)
    class(result) <- c("LRscreening", "list")
    result
}
#################################################################################
## S3 summary method for LRscreening objects produced by logisticScreening()
#################################################################################
#' Title in Title Case
#'
#' \code{summary.LRscreening} Short description.
#'
#' Longer description,
#'
#' @param Arguments.
#'
#' @return Contents of "Value"
#'
#' @author Steve Gutreuter, \email{sgutreuter@@gmail.com}
#'
#' @examples
#'
summary.LRscreening <- function(obj){
    if(!("LRscreening" %in% class(obj))) stop("obj not LRscreening class")
    cat("Call:\n")
    print(obj$Call)
    cat("\n\nLogistic regression model summary:")
    print(summary(obj$ModelFit))
    cat("\nOut-of-sample sensitivity and specificity\nscreening at p.hat >= p.threshold:\n\n")
    print(obj$CrossValPerf)
}
#################################################################################
## S3 print method for LRscreening objects produced by logisticScreening()
#################################################################################
#' Title in Title Case
#'
#' \code{print.LRScreening} Short description.
#'
#' Longer description,
#'
#' @param Arguments.
#'
#' @return Contents of "Value"
#'
#' @author Steve Gutreuter, \email{sgutreuter@@gmail.com}
#'
print.LRscreening <- function(obj){
    if(!("LRscreening" %in% class(obj))) stop("obj not LRscreening class")
    cat("Out-of-sample sensitivity and specificity\nscreening at p.hat >= p.threshold:\n\n")
    print(obj$CrossValPerf)
}
#################################################################################
## S3 plot method for LRscreening objects produced by logisticScreening()
#################################################################################
#' Title in Title Case
#'
#' \code(plot.LRscreening) Short description.
#'
#' A wrapper function for \code{lattice::xyplot} customized to produce a plot
#' of Receiver Operating Characteristic curves.
#'
#' @param Arguments.
#'
#' @return A \code{lattice} graphical object
#'
#' @author Steve Gutreuter, \email{sgutreuter@@gmail.com}
#'
#' @examples
#'
#' @export
plot.LRscreening <- function(obj, main = "Receiver Operating Characteristics",
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
