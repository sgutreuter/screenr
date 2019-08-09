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


#' Inverse of the Link for a Linear Predictor
#'
#' \code{inverseLink}
#'
#' Returns the inverse of logit, cloglog and probit link functions for a linear
#' predictor
#'
#' @param link Currently one of "logit", "cloglog" or "probit"
#' @param lp The linear predictor
#'
#' @return The inverse of the link function for the linear predictor.
#'
#' @author Steve Gutreuter, \email{sgutreuter@@gmail.com}
inverseLink <- function(link, lp){
    if(link == "logit"){
        p <- locfit::expit(lp)
    }else{
        if(link == "cloglog"){
            p <- 1 - exp(-exp(lp))
        } else {
            p <- pnorm(lp)
        }
    }
    p
}


#' #' An S3 Plot Method
#'
#' \code{plot.ROC} An S3 plotting method for Receiver Operating Characteristic
#' objects created by function plot.Bandason().
#'
#' @param x A Bandason class object.
#' @export
plot.ROC <- function(x){
    if(!("ROC" %in% class(x))) stop("x must be an ROC object")
    x$FPP <- 1 - x$specificity
    plot(x$FPP, x$sensitivity, type = "S", ylim = c(0,1),
         xlim = c(0, 1), ylab = "Sensitivity", xlab = "1 - Specificity",
         main = "Receiver Operating Characteristic")
    abline(a = 0, b = 1, lty = 2)
    text(x$sensitivity ~ x$FPP, labels = x$score, pos = 2)
}


#' Expected Number of Tests Required per Positive Result
#'
#' \code{testCounts} Expected tests per positive
#'
#' Compute the expected number of tests which need to be performed in order
#' to identify the first positive test result, and the expected number of
#' false positives among that number of tests.
#'
#' @param prev Proportion of the population expressing positive test results.
#' @param SensSpec A data frame containing columns '"sensitivity"' and
#' '"specificity"', or an object of class "Bandason" or class "binomscreenr".
#'
#' @return A data frame containing sensitivity, specificity, the expected
#' number of tests required to observe a single positive test result and,
#' among those, the expected number of false negatives per positive test
#' result.
#'
#' @author Steve Gutreuter, \email{sgutreuter@@gmail.com}
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
testCounts <- function(prev, SensSpec){
    if(!(prev > 0 & prev < 1)) stop("prev must be in (0,1)")
    if("binomscreenr" %in% class(SensSpec)){
        ss <- SensSpec[["CrossValPerf"]]
    } else {
        if("Bandason" %in% class(SensSpec)){
            ss <- SensSpec[["roc"]]
        } else {
            ss <- SensSpec
            if(!("sensitivity" %in% names(ss))) stop("No column 'sensitivity")
            if(!("specificity" %in% names(ss))) stop("No column 'specificity")
            if(any(ss[["sensitivity"]] < 0 | ss[["sensitivity"]] > 1))
                stop("sensitivity not in (0,1)")
            if(any(ss[["specificity"]] < 0 | ss[["specificity"]] > 1))
                stop("specificity not in (0,1)")
        }
    }
    prev <- prev[1]
    Etpp <- ((ss[["sensitivity"]] * prev) +
            (1 - ss[["specificity"]]) * (1 - prev)) / (ss[["sensitivity"]] * prev)
    Efn <- (1 - ss[["sensitivity"]]) * Etpp
    result <- data.frame(cbind(ss, Etpp, Efn))
    names(result)[ncol(ss) + c(1, 2)] <- c("E(Tests/Pos)",
                                           "E(FalseNegs/Pos)")
    result
}

################################   END of FILE   ################################
