################################################################################
##     R Script: confint-methods.R
##
##      Package: screenr
##
##  Description: S3 methods for confint of screenr objects
##
##       Author: Steve Gutreuter
##               sgutreuter@gmail.com
################################################################################


## Function confint.logreg_screenr
##
#' An S3 Method to Compute Confidence Limits from \code{logreg_screenr} Objects
#'
#' @description
#' \code{confint.logreg_screenr} returns the logistic model parameter estimates
#' and their profile-likelihood confidence limits from \code{logreg_screenr}-class
#' objects.
#'
#' @param object an object of class \code{logreg_screenr}.
#'
#' @param intercept (logical) retain (\code{TRUE}, default) or drop
#' (\code{FALSE}) the intercept coefficients.
#'
#' @param or return odds ratios if \verb{TRUE}. Default: FALSE (returns
#' logit-scale coefficients).
#'
#' @param conf_level confidence level for profile-likelihood confidence
#' intervals. Default: 0.95.
#'
#' @param digits number of decimal places to be printed. Default: 4.
#'
#' @param ... optional arguments passed to \code{predict} methods.
#'
#' @return \code{confint.logreg_screenr} returns a dataframe containing the
#' estimated coefficients (or odds ratios) and their profile-likelihood
#' lower and upper confidence limits (lcl and ucl, respectively).
#'
#' @examples
#' attach(uniobj2)
#' confint(uniobj2, or = TRUE)
#'
#' @importFrom stats confint
#' @export
confint.logreg_screenr <- function(object, ..., intercept =  TRUE, or = FALSE,
                                conf_level =  0.95, digits = 4){
    coef_ <- object[["ModelFit"]][["coefficients"]]
    plci_ <- confint(get_what(object, what = "ModelFit"), level = conf_level)
    if(intercept == FALSE){
        coef_ <- coef_[-1]
        plci_ <- plci_[-1, ]
        }
    if(or == TRUE ){
        coef_ <- exp(coef_)
        plci_ <- exp(plci_)
    }
    result <- data.frame(covariate = row.names(plci_),
                         estimate = round(coef_, digits = digits),
                         lcl = round(plci_[, 1], digits = digits),
                         ucl = round(plci_[, 2], digits = digits))
    row.names(result) <- NULL
    result
}


## Function confint.geeglm
##
#' An S3 Method to Compute Confidence Limits from \code{geepack::geeglm} Objects.
#'
#' @description \code{confint.geeglm} returns the logistic model coefficients
#' estimates and their  normal-theory Wald-type confidence
#' limits from objects produced by \code{geepack::geeglm}.
#'
#' @param object an object of class \code{geeglm}.
#'
#' @param conf_level confidence level in (0, 1).
#'
#' @param ... optional arguments passed to \code{predict} methods.
#'
#' @return \code{confint.geeglm} returns a dataframe containing the estimated
#' model coefficients and their lower and upper confidence limits, \code{lcl}
#' and \code{ucl}, respectively.
#'
#' @importFrom stats qnorm coef
#' @export
confint.geeglm <- function(object, ..., conf_level = 0.95) {
    x <- coef(summary(object))
    qn <- qnorm((1 + conf_level) / 2)
    result <- data.frame(covariate = rownames(x),
                         estimate = x$Estimate,
                         lcl = x$Estimate - qn * x$Std.err,
                         ucl = x$Estimate + qn * x$Std.err)
    result
}


## Function confint.gee_screenr
##
#' An S3 Method to Compute Confidence Limits from \code{gee_screenr} Objects
#'
#' @description
#' \code{confint.logreg_screenr} returns the logistic model parameter estimates
#' and their and Wald-type confidence limits from \code{gee_screenr}-class
#' objects.
#'
#' @param object an object of class \code{gee_screenr}.
#'
#' @param intercept (logical) retain (\code{TRUE}, default) or drop
#' (\code{FALSE}) the intercept coefficients.
#'
#' @param or return odds ratios if \verb{TRUE}. Default: FALSE (returns
#' logit-scale coefficients).
#'
#' @param conf_level confidence level for normal-theory Wald-type confidence
#' intervals. Default: 0.95.
#'
#' @param digits number of decimal places to be printed. Default: 4.
#'
#' @param ... optional arguments passed to \code{predict} methods.
#'
#' @return \code{confint.gee_screenr} returns a dataframe containing the
#' estimated coefficients (or odds ratios) and their Wald-type
#' lower and upper confidence limits (lcl and ucl, respectively).
#'
#' @examples
#' attach(uniobj3)
#' confint(uniobj3, or = TRUE)
#'
#' @importFrom stats confint
#' @export
confint.gee_screenr <- function(object, ..., intercept =  TRUE, or = FALSE,
                                conf_level =  0.95, digits = 4){
    x <- confint(get_what(object, what = "ModelFit"), conf_level = conf_level)
    if(intercept == FALSE){
        x <- x[-1, ]
        }
    if(or == TRUE ){
        x$estimate <- exp(x$estimate)
        x$lcl <- exp(x$lcl)
        x$ucl <- exp(x$ucl)
    }
    x$estimate <- round(x$estimate, digits = digits)
    x$lcl <- round(x$lcl, digits = digits)
    x$ucl <- round(x$ucl, digits = digits)
    x
}
