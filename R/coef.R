


## Function coef.glmpathScreenr
##
#' A method to extract the estimated coefficients from \code{glmpathScreenr} objects.
#'
#' @param object an object of class \code{glmpathScreenr}.
#'
#' @param intercept (logical) retain (\code{TRUE}, default) or drop
#' (\code{FALSE}) the intercept coefficients.
#'
#' @param or return odds ratios if \verb{TRUE}; logit-scale coefficients
#' are the default.
#'
#' @return a \emph{p} x 2 matrix containing the estimated coefficients from the AIC-
#' and BIC-best logistic regression models, where \emph{p} is the number of coefficients.
#'
#' @examples
#' attach(uniobj1)
#' coef(uniobj)
#'
#' @export
coef.glmpathScreenr <- function(object, intercept = TRUE, or = FALSE){
    if(!("glmpathScreenr" %in% class(object)))
        stop("object not glmpathScreenr class")
    coef_ <- data.frame(rbind(object$isResults$minAIC$Coefficients,
                              object$isResults$minBIC$Coefficients))
    rownames(coef_) <- c("AIC-best model", "BIC-best model")
    if (intercept == FALSE) coef_ <- coef_[, -1]
    if (or == TRUE) coef_ <- exp(coef_)
    coef_ <- t(as.matrix(coef_))
    coef_
}


## Function coef.logisticScreenr
##
#' \code{coef.logisticScreenr} is an S3 coef method for \code{logisticScreenr} objects.
#'
#' @param x an object of class \code{logisticScreenr}
#'
#' @param intercept (logical) retain (\code{TRUE}, default) or drop
#' (\code{FALSE}) the intercept coefficients.
#'
#' @param or return odds ratios if \verb{TRUE}; logit-scale coefficients
#' are the default.
#'
#' @details
#' Extracts the estimated coefficients from \code{logisticScreenr} objects.
#'
#' @return A numeric vector containing the estimated coefficients on the logit
#' scale.
#'
#' @examples
#' load(uniobj2)
#' class(uniobj2)
#' coef(uniobj2)
#'
#' @export
coef.logisticScreenr <- function(x, intercept =  TRUE, or = FALSE){
    stopifnot(class(x) == "logisticScreenr")
    coef_ <- x[["ModelFit"]][["coefficients"]]
    if(intercept == FALSE) coef_ <- coef_[-1]
    if(or == TRUE ) coef_ <- exp(coef_)
    coef_
}
