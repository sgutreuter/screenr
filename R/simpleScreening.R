#################################################################################
##       R PROGRAM: simpleScreening.R
##
##         PROJECT: R functions for HIV screening tool development
##
##      WRITTEN BY: Steve Gutreuter
##                  E-mail:  sgutreuter@gmail.com
##
#################################################################################

## Function simpleScreenr
##
#' \code{simpleScreenr} is used to produce overly simple test-screening tools.
#'
#' @param formula an object of class \code{\link[stats]{formula}} defining the
#' testing outcome and predictor covariates.
#' @param data the "training" sample; a data frame containing the testing outcome
#' and predictive covariates to be used for testing screening.  The testing
#' outcome must be binary (0,1) indicating negative and positive test results,
#' respectively, or logical (\verb{TRUE}/\verb{FALSE}), and the screening scores
#' are the row-wise sums of the values of those covariates.  The covariates are
#' typically binary (0 = no, 1 = yes) responses to questions, but the responses
#' may also be ordinal numeric values.
#'
#' @return An object of class \code{simpleScreenr} containing the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{Prevalence}}{Prevalence of the test condition in the training sample.}
#' \item{\code{ISroc}}{An object of class \code{\link[pROC]{roc}} containing
#' the "in-sample" (overly-optimistic) receiver operating characteristics,
#' and additional functions for use with this object are available in the
#' \code{pROC} package.}
#' \item{\code{Scores}}{The training sample, including the scores.}
#' }
#'
#' @details
#' \code{simpleScreenr} computes the in-sample (\emph{overly optimistic})
#' performances for development of a very simple test screening tool implementing
#' the method of Bandason et al. (2016).  \code{simpleScreener} is not optimized
#' and is intended only for comparision with \code{glmpathScreenr} or
#' \code{logisticScreenr}, either of which will almost certainly  out-perform
#' \code{simpleScreening}.
#'
#' @references
#' Bandason T, McHugh G, Dauya E, Mungofa S, Munyati SM, Weiss HA, Mujuru H,
#' Kranzer K, Ferrand RA. Validation of a screening tool to
#' identify older children living with HIV in primary care facilities in high
#' HIV prevalence settings. AIDS. 2016;30(5):779-785
#' \url{http://dx.doi.org/10.1097/QAD.0000000000000959}
#'
#' Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez J-C,
#' Müller M. \code{pROC}: An open-source package for \code{R} and S+ to
#' analyze and compare ROC curves. BMC Bioinformatics. 2011;12(77):1-8.
#' \url{http://doi.org/10.1186/1471-2105-12-77}
#'
#'
#' @examples
#' data(unicorns)
#' toosimple <- simpleScreenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5,
#'                           data = unicorns)
#' summary(toosimple)
#' \dontrun{testCounts(toosimple)}
#'
#' @seealso \code{\link{glmpathScreenr}}, \code{\link{logisticScreenr}}
#' @import pROC
#' @importFrom stats model.response complete.cases
#' @export
simpleScreenr <- function(formula, data){
    warning("WARNING! WARNING! WARNING! simpleScreenr is suboptimal and is provided only for comparison with other methods." )
    mf <- match.call(expand.dots = FALSE)
    call <- mf
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    dat <- data[complete.cases(data[, names(mf), ]), ]
    y <- model.response(mf, "numeric")
    if(!all(y %in% c(0, 1))) stop("Response variable must be binary (0, 1)")
    prev <- mean(y, na.rm = TRUE)
    preds <- mf[, -1]
    npreds <- dim(preds)[2]
    score <- apply(preds, 1, sum)
    is.roc <- pROC::roc(y, score, auc = TRUE, direction = "<")
    scores <- cbind(dat, score = score)
    result <- list(Call = call,
                   Prevalence = prev,
                   ISroc = is.roc,
                   Scores = scores)
    class(result) <- "simplescreenr"
    invisible(result)
}


## Function summary.simpleSreenr
##
#' \code{summary.simpleScreenr} is a summary method for \code{simpleScreenr} objects
#'
#' @param object an object of class \code{simpleScreenr}.
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @return Nothing.  Thresholds, specificities and sensitivities are printed as
#' a side effect.
#'
#' @seealso \code{\link{getROC}}
#' @export
summary.simpleScreenr <- function(object, ...){
    if(!("simpleScreenr" %in% class(object)))
        stop("object not a simpleScreenr object")
    cat("Call:\n")
    print(object$Call)
    cat("\nPrevalence (In-sample prevalence of condition):\n")
    print(object$Prevalence)
    cat("\nReceiver Operating Characteristics:\n")
    auc <- round(as.numeric(object$ISroc$auc), digits = 4)
    cat(paste("\nIn-sample area under the ROC curve: ",
              auc, "\n", sep = ""))
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
#' @param print_ci logical indicator for returning a dataframe of
#' numerical values (default = \verb{TRUE}).
#'
#' @param conf_level confidence level in the interval (0,1). Default is 0.95
#' producing 95\% confidence intervals.
#'
#' @param bootreps numeric-valued number of bootstrap replication for estimation
#' of 95\% confidence intervals.
#'
#' @param ... additional arguments for base \verb{\link{plot}} or passed to \verb{\link{plot.roc}} and friends.
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
#' toosimple <- simpleScreenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5,
#'                           data = unicorns)
#' plot(toosimple, type = "S")
#' @importFrom graphics plot
#' @export
plot.simpleScreenr <- function(x, plot_ci = TRUE, print_ci = TRUE,
                               conf_level = 0.95, bootreps = 2000,...){
    if(!class(x) == "simpleScreenr") stop("x is not a simpleScreenr object")
    plt <- plot(x$ISroc, print.auc = TRUE, ...)
    if(plot_ci | print_ci){
        ciplt <- pROC::ci.thresholds(x$ISroc, boot.n = bootreps,
                                     progress = "none",
                                     conf.level = conf_level,
                                     thresholds = "local maximas")
        }
    if(print_ci){
        threshold <- as.numeric(rownames(ciplt$sensitivity)) + 0.5
        citable <- data.frame(cbind(threshold, ciplt$sensitivity,
                                    ciplt$specificity))
        names(citable) <- c("threshold", "se.low", "se.median",
                            "se.high", "sp.low", "sp.median", "sp.high")
        row.names(citable)  <- 1:length(threshold)
    }
    if(plot_ci) plot(ciplt)
    if(print_ci) return(citable)
}


## Function print.simpleScreenr
##
#' \code{print.simpleScreenr} is an S3 print method for \code{simpleScreenr} objects.
#'
#' @param x an object of class \code{simplescreenr}.
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @return Nothing. Thresholds, specificities and sensitivities are printed as a
#' side effect.
#'
#' @details
#' Plots an ROC curve.
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
#' @seealso \code{\link{getROC}} and \code{\link{plot.simpleScreenr}}
#' @export
print.simpleScreenr <- function(x, ...){
    if(!("simpleScreenr" %in% class(x))) stop("x not a simpleScreenr object")
    cat("\nIn-sample (overly optimistic) sensitivity and specificity:\n")
    df_ <- pROC::coords(x$ISroc, transpose = FALSE)
    df_["threshold"] <- df_["threshold"] + 0.5
    print(df_)
}

################################   END of FILE   ################################
