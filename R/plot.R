################################################################################
##     R Script: plot.R
##
##      Package: screenr
##
##  Description: S3 methods for plotting
##
##       Author: Steve Gutreuter
##               sgutreuter@gmail.com
################################################################################


## Function plot.easyTool
##
#' \code{plot.easyTool} is a plotting method for \code{easyTool}
#' objects.
#' @param x an object of class \code{easyTool}.
#' @param plot_ci (logical) plot confidence intervals if \verb{TRUE}.
#' @param conf_level confidence level
#' @param bootreps the number of bootstrap replications for estimation of
#' confidence intervals.
#' @param ... any additional arguments passed to \code{pROC::plot.roc} or
#' \code{pROC::lines.roc}.
#'
#' @importFrom graphics legend plot
#'
#' @return
#'
#' This function produces a plot as a side effect and (optionally)
#' returns a dataframe containing sensitivities, specificities and their
#' lower and upper confidence limits for threshold values of Pr(response = 1).
#'
#' @details \code{plot.easyTool} is an enhanced convenience wrapper for
#' \code{`pROC::plot.roc`}.  The table is useful for identifying the
#' minimum predicted response probabilities associated with particular
#' sensitivities.  The sensitivities and specificities are the coordinates at
#' change points in the cross-validated ROC curve, and the values of Threshold_Pr
#' are the values of lower bound of the predicted probability that acheives those
#' sensitivities and specificities.  For example, if Threshold_P = 0.002, then
#' classifiction as positive of all those subjects for whom the predicted response
#' probability is greater than or equal to 0.002 will achieve the corresponding
#' sensitivity and specificity at the prescribed confidence level in new
#' data drawn from the same distribution as the training sample.
#'
#' @references
#' Fawcett T. An introduction to ROC analysis. Pattern Recognition Letters. 2006.
#' 27(8):861-874.
#' \url{https://doi.org/10.1016/j.patrec.2005.10.010}
#'
#' Linden A. Measuring diagnostic and predictive accuracy in disease
#' management: an introduction to receiver operating characteristic (ROC) analysis.
#' Journal of Evaluation in Clinical Practice. 2006; 12(2):132-139.
#' \url{https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1365-2753.2005.00598.x}
#'
#' Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez J-C, Muller M.
#' pROC: an open-source package for R and S+ to analyze and compare ROC curves.
#' BMC Bioinformatics 2011; 12:77. \url{https://www.biomedcentral.com/1471-2105/12/77}
#'
#' @examples
#' attach(uniobj1)
#' tool <- easyTool(uniobj1, max = 3, crossval = TRUE)
#' plot(tool)
#' @importFrom graphics legend plot lines
#' @import pROC
#' @export
plot.easyTool <- function(x, ..., plot_ci = TRUE,
                                 conf_level = 0.95, bootreps = 2000){
    if(!("easyTool" %in% class(x)))
            stop("Object not easyTool class")
    stopifnot(conf_level > 0 & conf_level < 1)
    roc_  <- x$ROC
    pROC::plot.roc(roc_, print.auc = TRUE, ci = FALSE, ...)
    if(plot_ci){
        ciplt <- pROC::ci.thresholds(roc_,
                                     boot.n = bootreps,
                                     progress = "text",
                                     conf.level = conf_level,
                                     thresholds = "local maximas")
    }
    if(plot_ci) plot(ciplt)
}


## Function plot.glmpathScreenr
##
#' \code{plot.glmpathScreenr} is a plotting method for \code{glmpathScreenr}
#' objects.
#'
#' @param x an object of class \code{glmpathScreenr}.
#' @param plot_ci (logical) plot confidence intervals if \verb{TRUE}.
#' @param model (character) select either the model which produced the
#' minimum AIC (\verb{"minAIC"}) or minimum BIC (\verb{"minBIC"}).
#' @param conf_level confidence level
#' @param bootreps the number of bootstrap replications for estimation of
#' confidence intervals.
#' @param ... any additional arguments passed to \code{pROC::plot.roc} or
#' \code{pROC::lines.roc}.
#'
#' @return This function produces a plot as a side effect.
#'
#' @details \code{plot.glmpathScreenr} is an enhanced convenience wrapper for
#' \code{`pROC::plot.roc`}.  The table is useful for identifying the
#' minimum predicted response probabilities associated with particular
#' sensitivities.  The sensitivities and specificities are the coordinates at
#' change points in the cross-validated ROC curve, and the values of Threshold_Pr
#' are the values of lower bound of the predicted probability that acheives those
#' sensitivities and specificities.  For example, if Threshold_P = 0.002, then
#' classifiction as positive of all those subjects for whom the predicted response
#' probability is greater than or equal to 0.002 will achieve the corresponding
#' sensitivity and specificity at the prescribed confidence level in new
#' data drawn from the same distribution as the training sample.
#'
#' @references
#' Fawcett T. An introduction to ROC analysis. Pattern Recognition Letters. 2006.
#' 27(8):861-874.
#' \url{https://doi.org/10.1016/j.patrec.2005.10.010}
#'
#' Linden A. Measuring diagnostic and predictive accuracy in disease
#' management: an introduction to receiver operating characteristic (ROC) analysis.
#' Journal of Evaluation in Clinical Practice. 2006; 12(2):132-139.
#' \url{https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1365-2753.2005.00598.x}
#'
#' Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez J-C, Muller M.
#' pROC: an open-source package for R and S+ to analyze and compare ROC curves.
#' BMC Bioinformatics 2011; 12:77. \url{https://www.biomedcentral.com/1471-2105/12/77}
#'
#' @examples
#' attach(uniobj1)
#' plot(uniobj1, model = "minAIC")
#' @importFrom graphics legend plot lines
#' @import pROC
#' @export
plot.glmpathScreenr <- function(x, ...,  plot_ci = TRUE, model = "minAIC",
                                conf_level = 0.95, bootreps = 2000){
    if(!("glmpathScreenr" %in% class(x)))
            stop("Object not glmpathScreenr class")
    stopifnot(conf_level > 0 & conf_level < 1)
    if(!model %in% c("minAIC", "minBIC"))
        stop("Specify 'minAIC' or 'minBIC' for model")
    cvROC <- x$cvResults[[model]][["ROC"]]
    isROC <- x$isResults[[model]][["ROC"]]

    pROC::plot.roc(cvROC, print.auc = TRUE, ci = FALSE, ...)
    if(plot_ci){
        ciplt <- pROC::ci.thresholds(cvROC,
                                     boot.n = bootreps,
                                     progress = "text",
                                     conf.level = conf_level,
                                     thresholds = "local maximas")
    }
    if(plot_ci) plot(ciplt)
    pROC::lines.roc(isROC, lty = 3)
    legend("bottomright", legend = c("cross-validated", "in-sample"),
           lty = c(1, 3), lwd = c(2, 2))
}


## Function plot.logisticScreenr
##
#' \code{plot.logisticScreenr} is an S3 plot method for \code{logisticScreenr} objects,
#'
#' @param x an object of class \code{logisticScreenr}.
#' @param plot_ci logical indicator for plotting point-wise confidence
#' intervals at the locally maximum subset of coordinates for
#' on sensitivity and specificity (default = \verb{TRUE}). See also
#' \code{\link[pROC]{ci.thresholds}}.
#' @param conf_level confidence level in the interval (0,1). Default is 0.95
#' producing 95\% confidence intervals
#' @param bootreps number of bootstrap replications for estimation of confidence
#' (default = 2000).
#' @param ... additional arguments passed to \code{\link[pROC]{plot.roc}} and friends.
#'
#' @return This function produces a plot as a side effect.
#'
#' @details
#' Plot cross-validated (out-of-sample) ROC curve with pointwise confidence
#' intevals along with the overly optimistic in-sample ROC curve.
#'
#' @references
#' Fawcett T. An introduction to ROC analysis. Pattern Recognition Letters. 2006.
#' 27(8):861-874.
#' \url{https://doi.org/10.1016/j.patrec.2005.10.010}
#'
#' Linden A. Measuring diagnostic and predictive accuracy in disease
#' management: an introduction to receiver operating characteristic (ROC) analysis.
#' Journal of Evaluation in Clinical Practice. 2006; 12(2):132-139.
#' \url{https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1365-2753.2005.00598.x}
#'
#' Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez J-C, Muller M.
#' pROC: an open-source package for R and S+ to analyze and compare ROC curves.
#' BMC Bioinformatics 2011; 12:77. \url{https://www.biomedcentral.com/1471-2105/12/77}
#'
#' @examples
#' attach(uniobj2)
#' plot(uniobj2)
#'
#' @importFrom graphics legend plot
#' @export
plot.logisticScreenr <- function(x, ..., plot_ci = TRUE, print = TRUE,
                              conf_level = 0.95, bootreps = 2000){
    if(!class(x) == "logisticScreenr") stop("x is not a logisticScreenr object")
    stopifnot(conf_level > 0 & conf_level < 1)
    plot(x$CVroc, print.auc = TRUE, ci = FALSE, ...)
    if(plot_ci | print){
        ciplt <- pROC::ci.thresholds(x$CVroc,
                                     boot.n = bootreps,
                                     progress = "text",
                                     conf.level = conf_level,
                                     thresholds = "local maximas")
    }
    if(plot_ci) plot(ciplt)
    pROC::lines.roc(x$ISroc, lty = 3)
    legend("bottomright", legend = c("cross-validated", "in-sample"),
           lty = c(1, 3), lwd = c(2, 2))
}


## Function plot.simpleScreenr
##
#' \code{plot.simpleScreenr} is a plot method for \code{simpleScreenr} objects.
#'
#' Plot ROC curve with pointwise 95% confidence
#' intevals on sensitivity and specificity and (optionally) returns a dataframe
#' containing numerical values.
#'
#' @param x an object of class \code{simpleScreenr}.
#'
#' @param plot_ci logical indicator for plotting point-wise confidence
#' intervals at the locally maximum subset of coordinates for
#' on sensitivity and specificity (default = \verb{TRUE}). See also
#' \code{\link[pROC]{ci.thresholds}}.
#'
#' @param conf_level confidence level in the interval (0,1). Default is 0.95
#' producing 95\% confidence intervals.
#'
#' @param bootreps numeric-valued number of bootstrap replication for estimation
#' of 95\% confidence intervals.
#'
#' @param ... additional arguments for \verb{\link{plot}} or passed to \verb{\link{plot.roc}} and friends.
#'
#' @return This function produces a plot as a side effect, and (optionally)
#' returns a dataframe dataframe containing medians and
#' bootstrap confidence limits of sensitivity and specificity.
#'
#' @references
#' Fawcett T. An introduction to ROC analysis. Pattern Recognition Letters. 2006.
#' 27(8):861-874.
#' \url{https://doi.org/10.1016/j.patrec.2005.10.010}
#'
#' Linden A. Measuring diagnostic and predictive accuracy in disease
#' management: an introduction to receiver operating characteristic (ROC) analysis.
#' Journal of Evaluation in Clinical Practice. 2006; 12(2):132-139.
#' \url{https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1365-2753.2005.00598.x}
#'
#' Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez J-C, Muller M.
#' pROC: an open-source package for R and S+ to analyze and compare ROC curves.
#' BMC Bioinformatics 2011; 12:77. \url{https://www.biomedcentral.com/1471-2105/12/77}
#'
#' @examples
#' data(unicorns)
#' too_simple <- simpleScreenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
#'                           data = unicorns)
#' plot(too_simple)
#' @importFrom graphics plot
#' @export
plot.simpleScreenr <- function(x, ..., plot_ci = TRUE, conf_level = 0.95,
                               bootreps = 2000){
    if(!class(x) == "simpleScreenr") stop("x is not a simpleScreenr object")
    plt <- plot(x$ISroc, print.auc = TRUE, ...)
    if(plot_ci){
        ciplt <- pROC::ci.thresholds(x$ISroc, boot.n = bootreps,
                                     progress = "none",
                                     conf.level = conf_level,
                                     thresholds = "local maximas")
        }
    if(plot_ci) plot(ciplt)
}
