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
