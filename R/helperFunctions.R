#################################################################################
##       R PROGRAM: helperFunctions.R
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
#' Estimate the required number of positive and negative test results for
#' hypothesis testing against the joint null hypothesis H0: Sens = SnsCrit and
#' Spec = SpcCrit.
#'
#' @param Sens A vector of anticipated sensitivities (0 < Sens < 1).
#' @param Spec A vector of anticipated specificities (0 < Spec < 1).
#' @param SnsCrit Critical value for sensitivity.
#' @param SpcCrit Critical value for specificity.
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
#' nSensSpec(sens, spec, SnsCrit = 0.90, SpcCrit = 0.70)
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
#' Computes sensitivity and specificity of a test.
#' @param x A 2 x 2 table, with the numbers of negative test results appearing
#' first in both rows and columns.
#' @return A list containing components sensitivity and specificity.
#' @examples
#' TestResults <- ordered(c(0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0))
#' TrueResults <- ordered(c(0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0))
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


#' Inverses of Binomial Regression Link Functions
#'
#' Returns the inverse of logit, cloglog and probit link functions for a linear
#' predictor
#'
#' @param link Character link function (one of "logit", "cloglog" or "probit")
#' @param lp Numeric linear predictor
#'
#' @return The inverse of the link function for the linear predictor.
#' @export
inverseLink <- function(link, lp){
    if(!link %in% c("logit", "cloglog", "probit")) stop("Bad link specification")
    if(link == "logit"){
        p <- exp(lp) / (1 + exp(lp))
    }else{
        if(link == "cloglog"){
            p <- 1 - exp(-exp(lp))
        } else {
            p <- pnorm(lp)
        }
    }
    p
}


#' Plot Receiver Operating Characteristics
#'
#' @param x an object of class \code{binomscreener} \code{simplescreenr}.
#' @param main plot title.
#' @param ... arguments to be passed to or from \code{lattice::xyplot}.
#' @return A \code{lattice} graphical object.
#' @export
plotROC <- function(x, main = "Receiver Operating Characteristic", ...){
    if(!((class(x) == "binomscreenr") | (class(x) == "simplescreenr")))
        stop("x not a 'binomscreenr' or 'simplescreenr' object")
    if(class(x) == "binomscreenr"){
            nrows <- dim(x$CrossValPerf)[1]
            d.lty <- c(2, 1)
            dfrm <- data.frame(p.threshold =
                                   rep(x$CrossValPerf$p.threshold, 2),
                               grp = c(rep("In-sample", nrows),
                                       rep("Out-of-sample", nrows)),
                               sensitivity  = c(x$InSamplePerf$sensitivity,
                                                x$CrossValPerf$sensitivity),
                               FPP = 1 - c(x$InSamplePerf$specificity,
                                           x$CrossValPerf$specificity))
            d.key <- list(corner = c(1, 0), x = 0.98, y = 0.04,
                          text = list(c("In-sample (overly-optimistic)",
                                        "Out-of-sample")),
                          lines = list(type = c("l", "l"), col = rep("black", 2),
                                       lty = d.lty))
            res <- lattice::xyplot(sensitivity ~ FPP,
                                   group = grp,
                                   data = dfrm,
                                   panel = function(x, y, ...){
                                       panel.xyplot(x, y, ...)
                                       panel.abline(a = 0, b = 1)
                                   },
                                   main = main,
                                   type = rep("S", 2),
                                   col = rep("black", 2),
                                   lty = d.lty,
                                   ylab = "Sensitivity (%)",
                                   xlab = "1 - Specificity (%)",
                                   xlim = c(-0.02, 1),
                                   ylim = c(0, 1),
                                   key = d.key)
    } else {
        dfrm <- x$InSamplePerf
        dfrm$FPP <- 1 - dfrm$specificity
        res <- lattice::xyplot(sensitivity ~ FPP,
                               data = dfrm,
                               panel = function(x, y, ...){
                                   panel.xyplot(x, y, ...)
                                   panel.abline(a = 0, b = 1)
                                   panel.text(x, y,
                                              labels = as.character(data$score),
                                              pos = 2)
                               },
                               type = "S",
                               col = "black",
                               ylim = c(0, 1),
                               xlim = c(-0.02, 1),
                               ylab = "Sensitivity",
                               xlab = "1 - Specificity",
                               main = main)
    }
    res
}


#' Expected Number of Tests Required per Positive Test Result
#'
#' Compute the expected number of tests which need to be performed in order
#' to identify the first positive test result, and the expected number of
#' false positives among that number of tests.
#'
#' @param SensSpec A data frame containing columns "sensitivity" and
#' '"specificity"', or an object of class 'simplescreenr' or 'binomscreenr'.
#' @param prev Numeric proportion of the population expressing positive test
#' results.  \code{prev} is optional for class 'simplescreenr' and
#' 'binomscreenr' objects, and defaults to prevalence in the training sample
#' if not specified.
#'
#' @return A data frame containing sensitivity, specificity, the expected
#' number of tests required to observe a single positive test result and,
#' among those, the expected number of false negatives per positive test
#' result.
#'
#' @examples
#' data(unicorns)
#' unitool <- binomialScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5,
#'                              data = unicorns, Nfolds = 20,
#'                              p.threshold = c(seq(0.01, 0.10, by = 0.015),
#'                                              seq(.15, 0.95, by = 0.05)))
#' testCounts(unitool)
#'
#' @export
testCounts <- function(SensSpec = NULL, prev = NULL){
    if("binomscreenr" %in% class(SensSpec)){
        ss <- SensSpec[["CrossValPerf"]]
        if(is.null(prev)) prev <- SensSpec$Prevalence
    } else {
        if("simplescreenr" %in% class(SensSpec)){
            ss <- SensSpec[["InSamplePerf"]]
            if(is.null(prev)) prev <- SensSpec$Prevalence
        } else {
            ss <- SensSpec
            if(!("sensitivity" %in% names(ss))) stop("No column 'sensitivity")
            if(!("specificity" %in% names(ss))) stop("No column 'specificity")
            if(any(ss[["sensitivity"]] < 0 | ss[["sensitivity"]] > 1))
                stop("sensitivity not in (0,1)")
            if(any(ss[["specificity"]] < 0 | ss[["specificity"]] > 1))
                stop("specificity not in (0,1)")
            if(is.null(prev)) stop("Argument prev missing")
            if(!(prev > 0 & prev < 1)) stop("prev must be in (0,1)")
        }
    }
    Etpp <- ((ss[["sensitivity"]] * prev) +
            (1 - ss[["specificity"]]) * (1 - prev)) / (ss[["sensitivity"]] * prev)
    Efn <- (1 - ss[["sensitivity"]]) * Etpp
    result <- data.frame(cbind(ss, Etpp, Efn))
    names(result)[ncol(ss) + c(1, 2)] <- c("E(Tests/Pos)",
                                           "E(FalseNegs/Pos)")
    result
}

################################   END of FILE   ################################
