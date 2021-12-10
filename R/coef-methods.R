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
#' \code{logreg_screenr} objects.
#'
#' @return a \emph{p} x 2 matrix containing the estimated coefficients from the AIC-
#' and BIC-best logistic regression models, where \emph{p} is the number of coefficients.
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
    coef_
}


## Function coef.logreg_screenr
##
#' An S3 Method to Extract Coefficients from \code{logreg_screenr} Objects
#'
#' @param object an object of class \code{logreg_screenr}.
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
#' \code{coef.logreg_screenr} extracts the estimated coefficients from
#' \code{logreg_screenr} objects.
#'
#' @return A numeric vector containing the estimated coefficients on the logit
#' scale.
#'
#' @examples
#' attach(uniobj2)
#' class(uniobj2)
#' coef(uniobj2)
#'
#' @importFrom stats coef
#' @export
coef.logreg_screenr <- function(object, ..., intercept =  TRUE, or = FALSE){
    stopifnot(class(object) == "logreg_screenr")
    coef_ <- object[["ModelFit"]][["coefficients"]]
    if(intercept == FALSE) coef_ <- coef_[-1]
    if(or == TRUE ) coef_ <- exp(coef_)
    coef_
}
