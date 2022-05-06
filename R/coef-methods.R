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
#' \code{lasso_screenr} objects.  Regularization does not support estimation
#' of confidence limits.
#'
#' @return \code{coef.lasso_screenr} returns a dataframe containing the
#' estimated coefficients (or odds ratios) from the AIC- and BIC-best
#' logistic regression models, where \emph{p} is the number of coefficients.
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
#' from and profile-likelihood confidence limits from \code{logreg_screenr}-class
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
#' @details
#' \code{coef.logreg_screenr} extracts the estimated coefficients from
#' \code{logreg_screenr} objects.
#'
#' @return \code{coef.logreg_screenr} returns a dataframe containing the
#' estimated coefficients (or odds ratios) and their profile-likelihood
#' lower and upper confidence limits (lcl and ucl, respectively).
#'
#' @examples
#' attach(uniobj2)
#' coef(uniobj2, or = TRUE)
#'
#' @importFrom stats coef confint
#' @export
coef.logreg_screenr <- function(object, ..., intercept =  TRUE, or = FALSE,
                                conf_level =  0.95, digits = 4){
    stopifnot(class(object) == "logreg_screenr")
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
    coef_ <- data.frame(covariate = row.names(plci_),
                        estimate = round(coef_, digits = digits),
                        lcl = round(plci_[, 1], digits = digits),
                        ucl = round(plci_[, 2], digits = digits))
    row.names(coef_) <- NULL
    coef_
}
