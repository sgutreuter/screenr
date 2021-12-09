################################################################################
##     R Script: predict.R
##
##      Package: screenr
##
##  Description: S3 methods for prediction
##
##       Author: Steve Gutreuter
##               sgutreuter@gmail.com
################################################################################


## Function predict.glmpathScreenr
##
#' \code{predict.glmpathScreenr} is a prediction method for objects of class \code{glmpathScreenr}
#'
#' @param object an object of class \code{glmpathScreenr} produced by
#' \code{`glmpathScreenr`}.
#'
#' @param newdata new dataframe from which predicted probabilities of positive test results are desired.
#' The dataframe must contain values of the same response variables and covariates that were used
#' to obtain \code{obj}.
#'
#' @param ... optional arguments to \code{predict} methods.
#'
#' @return \code{predict.glmpathScreenr} returns (invisibly) a dataframe augmenting the complete cases
#' in \code{newdata} with the predicted probabilities of positive test results \code{phat_minAIC} and
#' \code{phat_minBIC} from the models that produced the minimum AIC and BIC, respectively.
#'
#' @details This method is a convenience wrapper for \code{`glmpath::predict.glmpath`}.
#'
#' @import glmpath
#' @importFrom stats predict
#' @examples
#' attach(uniobj1)
#' ## Get some new observations
#' new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
#'                         testresult = c(NA, NA),
#'                         Q1 = c(0, 0), Q2 = c(0, 0), Q3 = c(0, 1), Q4 = c(0, 0),
#'                         Q5 = c(0, 1), Q6 = c(0, 1), Q7 = c(0, 1))
#' ## Predict the probabilities of testing positive for the new subjects
#' predict(uniobj1, newdata = new_corns)
#' @export
predict.glmpathScreenr <- function(object =  NULL, ..., newdata = NULL){
    if(!is.data.frame(newdata)) stop("Specify a dataframe")
    if(!("glmpathScreenr" %in% class(object))) stop("object not a glmpathScreenr object")
    form  <- object$formula
    rname <- as.character(form[[2]])
    nd <- newdata
    nd[is.na(nd[[rname]]), rname] <- 0
    mf <- model.frame(form, nd)
    y <- mf[, 1]
    x <- as.matrix(mf[, -1])
    obj <- object$glmpathObj
    sAIC <- which(obj$aic == min(obj$aic))
    pAIC <- glmpath::predict.glmpath(obj, newx =  x,  newy = y,  s = sAIC, type = "response")
    sBIC <- which(obj$bic == min(obj$bic))
    pBIC <- glmpath::predict.glmpath(obj, newx =  x,  newy = y,  s = sBIC, type = "response")
    dat <- data.frame(cbind(y, x))
    names(dat)[1] <- rname
    dat[[rname]] <- newdata[[rname]]
    res <- data.frame(pAIC, pBIC)
    names(res) <- c("phat_minAIC", "phat_minBIC")
    res <- cbind(dat, res)
    invisible(res)
}


## Function predict.logisticScreenr
##
#' \code{predict.logisticScreenr} is an S3 prediction method for objects of class \code{logisticScreenr}
#'
#' @param object an object of class \code{logisticScreenr} produced by
#' \code{`logisticScreenr`}.
#'
#' @param newdata new dataframe from which predicted probabilities of positive test results are
#' desired. The dataframe must contain values of the same response variables and covariates that
#' were used to obtain \code{object}.
#'
#' @param ... optional arguments to \code{predict} methods.
#'
#' @return \code{predict.logisticScreenr} returns (invisibly) a dataframe augmenting
#' \code{newdata} with the predicted probabilities of positive test results \code{phat}.
#'
#' @details This method is a convenience wrapper for \code{`stats::predict.glm`}.
#'
#' @examples
#' attach(uniobj2)
#' ## Get some new observations
#' new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
#'                         testresult = c(NA, NA), Q1 = c(0, 0), Q2 = c(0, 0),
#'                         Q3 = c(0, 0), Q4 = c(0, 0), Q5 = c(0, 1), Q6 = c(0, 1 ),
#'                         Q7 = c(0, 1))
#' ## Predict the probabilities of testing positive for the new subjects
#' predict(uniobj2, newdata = new_corns)
#' @importFrom stats predict
#' @export
#'
#'
predict.logisticScreenr <- function(object = NULL, ...,  newdata = NULL){
    if(!is.data.frame(newdata)) stop("Specify a dataframe")
    if(!("logisticScreenr" %in% class(object))) stop("object not a logisticScreenr object")
    form <- object$formula
    rname <- as.character(form[[2]])
    nd <- newdata
    nd[is.na(nd[[rname]]), rname] <- 0
    res <- stats::predict.glm(object$ModelFit, newdata = nd, type = "response")
    res <- data.frame(newdata, data.frame(phat = res))
    invisible(res)
}
