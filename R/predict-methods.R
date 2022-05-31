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


## Function predict.lasso_screenr
##
#' A Prediction Method for \code{lasso_screenr}-Class Objects
#'
#' @description \code{predict.lasso_screenr} computes predicted probabilities of
#' positive test results from new data.
#'
#' @param object an object of class \code{lasso_screenr} produced by
#' \code{`lasso_screenr`}.
#'
#' @param newdata new dataframe from which predicted probabilities of positive
#' test results are desired.  The dataframe must contain values of the same
#' response variables and covariates that were used to obtain \code{obj}.
#'
#' @param ... optional arguments to \code{predict} methods.
#'
#' @return \code{predict.lasso_screenr} returns (invisibly) a dataframe
#' augmenting the complete cases in \code{newdata} with the predicted
#' probabilities of positive test results \code{phat_minAIC} and
#' \code{phat_minBIC} from the models that produced the minimum AIC and BIC,
#' respectively.
#'
#' @details This method is a convenience wrapper for
#' \code{`glmpath::predict.glmpath`}.
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
#' print(predict(uniobj1, newdata = new_corns))
#' @export
predict.lasso_screenr <- function(object =  NULL, ..., newdata = NULL){
    if(!is.data.frame(newdata)) stop("Specify a dataframe")
    if(!("lasso_screenr" %in% class(object)))
        stop("object not a lasso_screenr object")
    form  <- object$formula
    rname <- as.character(form[[2]])
    nd <- newdata
    nd[is.na(nd[[rname]]), rname] <- 0
    mf <- model.frame(form, nd)
    y <- mf[, 1]
    x <- as.matrix(mf[, -1])
    obj <- object$glmpathObj
    sAIC <- which(obj$aic == min(obj$aic))
    pAIC <- glmpath::predict.glmpath(obj, newx =  x,  newy = y,  s = sAIC,
                                     type = "response")
    sBIC <- which(obj$bic == min(obj$bic))
    pBIC <- glmpath::predict.glmpath(obj, newx =  x,  newy = y,  s = sBIC,
                                     type = "response")
    res <- data.frame(pAIC, pBIC)
    names(res) <- c("phat_minAIC", "phat_minBIC")
    res <- cbind(newdata, res)
    invisible(res)
}


## Function predict.logreg_screenr
##
#' A Prediction Method for \code{logreg_screenr}-Class Objects
#'
#' @description \code{predict.logreg_screenr} computes predicted probabilities
#' of positive test results from new data.
#'
#' @param object an object of class \code{logreg_screenr} produced by
#' \code{`logreg_screenr`}.
#'
#' @param newdata new dataframe from which predicted probabilities of positive
#' test results are desired. The dataframe must contain values of the same
#' response variables and covariates that were used to obtain \code{object}.
#'
#' @param ... optional arguments to \code{predict} methods.
#'
#' @return \code{predict.logreg_screenr} returns (invisibly) a dataframe
#' augmenting \code{newdata} with the predicted probabilities of positive test
#' results \code{phat}.
#'
#' @details This method is a convenience wrapper for \code{`stats::predict.glm`}.
#'
#' @examples
#' attach(uniobj2)
#' ## Get some new observations
#' new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
#'                         testresult = c(NA, NA), Q1 = c(0, 0), Q2 = c(0, 0),
#'                         Q3 = c(0, 0), Q5 = c(0, 1), Q6 = c(0, 1 ),
#'                         Q7 = c(0, 1))
#' ## Predict the probabilities of testing positive for the new subjects
#' print(predict(uniobj2, newdata = new_corns))
#' @importFrom stats predict
#' @export
predict.logreg_screenr <- function(object = NULL, ...,  newdata = NULL){
    if(!is.data.frame(newdata)) stop("Specify a dataframe")
    if(!("logreg_screenr" %in% class(object)))
        stop("object not a logreg_screenr object")
    form <- object$formula
    rname <- as.character(form[[2]])
    nd <- newdata
    nd[is.na(nd[[rname]]), rname] <- 0
    res <- stats::predict.glm(object$ModelFit, newdata = nd, type = "response")
    res <- data.frame(newdata, data.frame(phat = res))
    invisible(res)
}


## Function predict.easy_tool
##
#' A Prediction Method for \code{easy_tool}-Class Objects
#'
#' @description \code{predict.easy_tool} computes predicted simplified
#' screening scores from new data.
#'
#' @param object an object of class \code{easy_tool} produced by
#' \code{`easy_tool`}.
#'
#' @param newdata new dataframe from which predicted simplified screening
#' scores are desired. The dataframe must contain values of the same
#' response variables and covariates that were used to obtain \code{object}.
#'
#' @param ... optional arguments to \code{predict} methods.
#'
#' @return \code{predict.easy_tool} returns (invisibly) a dataframe
#' augmenting \code{newdata} with the predicted simplified test screening
#' scores \code{score}.
#'
#' @examples
#' new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
#'                         testresult = c(NA, NA), Q1 = c(0, 0), Q2 = c(0, 0),
#'                         Q3 = c(0, 0), Q4 = c(0, 0), Q5 = c(0, 1), Q6 = c(0, 1 ),
#'                         Q7 = c(0, 1))
#' attach(uniobj1)
#' et <- easy_tool(uniobj1, max = 3)
#' print(predict(et, newdata = new_corns))
#' @export
predict.easy_tool <- function(object =  NULL, ..., newdata =  NULL){
    if(!is.data.frame(newdata)) stop("Specify a dataframe")
    if(!("easy_tool" %in% class(object)))
        stop("object not an easy_tool object")
    xvarn <- object$varnames[-1]
    x = as.matrix(newdata[, xvarn])
    pred = x %*% object$QuestionWeights
    res = data.frame(cbind(newdata, score = pred))
    invisible(res)
}
