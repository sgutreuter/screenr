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

## Function plot.easy_tool
##
#' Plot ROC Curves from \code{easy_tool}-Class Objects
#'
#' @description
#' \code{plot.easy_tool} plots the \emph{k}-fold cross-validated
#' receiver-operating characteristics, including confidence intervals on the
#' combinations of the local maxima of sensitivity and specificity.
#'
#' @param x an object of class \code{easy_tool}.
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
#' This function produces a plot as a side effect and (optionally)
#' returns a dataframe containing sensitivities, specificities and their
#' lower and upper confidence limits for threshold values of Pr(response = 1).
#'
#' @details \code{plot.easy_tool} is an enhanced convenience wrapper for
#' \code{pROC::plot.roc}.
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
#' tool <- easy_tool(uniobj1, max = 3, crossval = TRUE)
#' plot(tool)
#' @importFrom graphics legend plot lines
#' @import pROC
#' @export
plot.easy_tool <- function(x, ..., plot_ci = TRUE,
                                 conf_level = 0.95, bootreps = 2000){
    if(!("easy_tool" %in% class(x)))
            stop("Object not easy_tool class")
    stopifnot(conf_level > 0 & conf_level < 1)
    roc_  <- x$ROC
    pROC::plot.roc(roc_, print.auc = TRUE, ci = FALSE, ...)
    if(plot_ci){
        ciplt <- pROC::ci.thresholds(roc_,
                                     boot.n = bootreps,
                                     conf.level = conf_level)
    }
    if(plot_ci) plot(ciplt)
}


## Function plot.lasso_screenr
##
#' Plot ROC Curves from \code{lasso_screenr}-Class Objects
#'
#' @description
#' \code{plot.lasso_screenr} plots the \emph{k}-fold cross-validated
#' receiver-operating characteristic, including confidence intervals on the
#' combinations of the local maxima of sensitivity and specificity.
#'
#' @param x an object of class \code{lasso_screenr}.
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
#' @details
#' Plot cross-validated (out-of-sample) ROC curve with pointwise confidence
#' intevals along with the overly optimistic in-sample ROC curve.
#' \code{plot.lasso_screenr} is an enhanced convenience wrapper for
#' \code{pROC::plot.roc}.
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
plot.lasso_screenr <- function(x, ...,  plot_ci = TRUE, model = c("minAIC", "minBIC"),
                                conf_level = 0.95, bootreps = 2000){
    if(!("lasso_screenr" %in% class(x)))
            stop("Object not lasso_screenr class")
    stopifnot(conf_level > 0 & conf_level < 1)
    model <- match.arg(model)
    cvROC <- x$cvResults[[model]][["ROC"]]
    isROC <- x$isResults[[model]][["ROC"]]

    pROC::plot.roc(cvROC, print.auc = TRUE, ci = FALSE, ...)
    if(plot_ci){
        ciplt <- pROC::ci.thresholds(cvROC,
                                     boot.n = bootreps,
                                     conf.level = conf_level)
    }
    if(plot_ci) plot(ciplt)
    pROC::lines.roc(isROC, lty = 3)
    legend("bottomright", legend = c("cross-validated", "in-sample"),
           lty = c(1, 3), lwd = c(2, 2))
}


## Function plot.logreg_screenr
##
#' Plot ROC Curves from \code{logreg_screenr}-Class Objects
#'
#' @description
#' \code{plot.logreg_screenr} plots the \emph{k}-fold cross-validated
#' receiver-operating characteristic, including confidence intervals on the
#' combinations of the local maxima of sensitivity and specificity.
#'
#' @param x an object of class \code{logreg_screenr}.
#' @param plot_ci logical indicator for plotting point-wise confidence
#' intervals at the locally maximum subset of coordinates for
#' on sensitivity and specificity (default = \verb{TRUE}). See also
#' \code{\link[pROC]{ci.thresholds}}.
#' @param conf_level confidence level in the interval (0,1). Default is 0.95
#' producing 95\% confidence intervals.
#' @param bootreps number of bootstrap replications for estimation of confidence
#' (default = 2000).
#' @param ... additional arguments passed to \code{\link[pROC]{plot.roc}} and friends.
#'
#' @return This function produces a plot as a side effect.
#'
#' @details
#' Plot cross-validated (out-of-sample) ROC curve with pointwise confidence
#' intevals along with the overly optimistic in-sample ROC curve.
#' \code{plot.lasso_screenr} is an enhanced convenience wrapper for
#' \code{pROC::plot.roc}.
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
plot.logreg_screenr <- function(x, ..., plot_ci = TRUE, conf_level = 0.95,
                                bootreps = 2000){
    if(!class(x) == "logreg_screenr") stop("x is not a logreg_screenr object")
    stopifnot(conf_level > 0 & conf_level < 1)
    plot(x$CVroc, print.auc = TRUE, ci = FALSE, ...)
    if(plot_ci){
        ciplt <- pROC::ci.thresholds(x$CVroc,
                                     boot.n = bootreps,
                                     conf.level = conf_level)
    }
    if(plot_ci) plot(ciplt)
    pROC::lines.roc(x$ISroc, lty = 3)
    legend("bottomright", legend = c("cross-validated", "in-sample"),
           lty = c(1, 3), lwd = c(2, 2))
}


## Function plot.simple_screenr
##
#' Plot ROC Curves from \code{simple_screenr}-Class Objects
#'
#' @description
#' \code{plot.simple_screenr} plots the \emph{k}-fold cross-validated
#' receiver-operating characteristic, including confidence intervals on the
#' combinations of the local maxima of sensitivity and specificity.
#'
#' Plot ROC curve with pointwise 95% confidence
#' intevals on sensitivity and specificity and (optionally) returns a dataframe
#' containing numerical values.
#'
#' @param x an object of class \code{simple_screenr}.
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
#' too_simple <- simple_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
#'                           data = unicorns)
#' plot(too_simple)
#' @importFrom graphics plot
#' @export
plot.simple_screenr <- function(x, ..., plot_ci = TRUE, conf_level = 0.95,
                               bootreps = 2000){
    if(!class(x) == "simple_screenr") stop("x is not a simple_screenr object")
    plt <- plot(x$ISroc, print.auc = TRUE, ...)
    if(plot_ci){
        ciplt <- pROC::ci.thresholds(x$ISroc, boot.n = bootreps,
                                     progress = "none",
                                     conf.level = conf_level)
        }
    if(plot_ci) plot(ciplt)
}
