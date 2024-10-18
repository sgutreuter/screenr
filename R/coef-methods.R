################################################################################
##     R Script: coef.R
##
##      Package: screenr
##
##  Description: S3 methods for coef of screenr objects
##
##       Author: Steve Gutreuter
##               sgutreuter@gmail.com
################################################################################


## Function coef.lasso_screenr
##
#' An S3 Method to Extract Coefficients from \code{lasso_screenr} Objects
#'
#' @description
#' \code{coef.lasso_screenr} returns the regularized logistic model parameter
#' estimates from the AIC- and BIC-best fits from \code{lasso_screenr}-class
#' objects.
#'
#' @param object an object of class \code{lasso_screenr}.
#'
#' @param intercept (logical) retain (\code{TRUE}, default) or drop
#' (\code{FALSE}) the intercept coefficients.
#'
#' @param or return odds ratios if \verb{TRUE}; logit-scale coefficients
#' are the default.
#'
#' @param ... optional arguments passed to \code{predict} methods.
#'
#' @details
#' \code{coef.lasso_screenr} extracts the estimated coefficients from
#' \code{lasso_screenr} objects.
#'
#' @examples
#' attach(uniobj1)
#' coef(uniobj1)
#' @importFrom stats coef
#' @export
coef.lasso_screenr <- function(object, ..., intercept = TRUE, or = FALSE){
    if(!("lasso_screenr" %in% class(object)))
        stop("object not lasso_screenr class")
    coef_ <- data.frame(rbind(object$isResults$minAIC$Coefficients,
                              object$isResults$minBIC$Coefficients))
    rownames(coef_) <- c("AIC-best model", "BIC-best model")
    if (intercept == FALSE) coef_ <- coef_[, -1]
    if (or == TRUE) coef_ <- exp(coef_)
    coef_ <- t(as.matrix(coef_))
    data.frame(coef_)
}


## Function coef.logreg_screenr
##
#' An S3 Method to Extract Coefficients from \code{logreg_screenr} Objects
#'
#' @description
#' \code{coef.logreg_screenr} returns the logistic model parameter estimates
#' from \code{logreg_screenr}-class objects.
#'
#' @param object an object of class \code{logreg_screenr}.
#'
#' @param intercept (logical) retain (\code{TRUE}, default) or drop
#' (\code{FALSE}) the intercept coefficients.
#'
#' @param or return odds ratios if \verb{TRUE}. Default: FALSE (returns
#' logit-scale coefficients).
#'
#' @param digits number of decimal places to be printed. Default: 4.
#'
#' @param ... optional arguments passed to \code{predict} methods.
#'
#' @seealso \code{\link[screenr]{confint.logreg_screenr}}
#'
#' @return \code{coef.logreg_screenr} returns a dataframe containing the
#' estimated coefficients (or odds ratios).
#'
#' @examples
#' attach(uniobj2)
#' coef(uniobj2, or = TRUE)
#'
#' @importFrom stats coef
#' @export
coef.logreg_screenr <- function(object, ..., intercept =  TRUE, or = FALSE,
                                digits = 4){
    if(!("logreg_screenr" %in% class(object)))
        stop("object not lasso_screenr class")
    coef_ <- object[["ModelFit"]][["coefficients"]]
    if(intercept == FALSE) coef_ <- coef_[-1]
    if(or == TRUE ) coef_ <- exp(coef_)
    coef_ <- data.frame(covariate = names(coef_),
                        estimate = round(coef_, digits = digits))
    row.names(coef_) <- NULL
    coef_
}
