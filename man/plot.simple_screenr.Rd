% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-methods.R
\name{plot.simple_screenr}
\alias{plot.simple_screenr}
\title{An S3 Method to Plot ROC Curves}
\usage{
\method{plot}{simple_screenr}(
  x,
  ...,
  plot_ci = TRUE,
  conf_level = 0.95,
  bootreps = 4000,
  print_auc = TRUE,
  partial_auc = c(0.8, 1),
  partial_auc_focus = c("sensitivity", "specificity"),
  partial_auc_correct = TRUE,
  type = c("l", "S")
)
}
\arguments{
\item{x}{an object of class \code{simple_screenr}.}

\item{...}{additional arguments for \verb{\link{plot}} or passed to
\verb{\link{plot.roc}} and friends.}

\item{plot_ci}{logical indicator for plotting point-wise confidence
intervals at the locally maximum subset of coordinates for
on sensitivity and specificity.  Default: \verb{TRUE}. See also
\code{\link[pROC]{ci.thresholds}}.}

\item{conf_level}{confidence level in the interval (0,1). Default is 0.95
producing 95\% confidence intervals. Default: \verb{TRUE}.}

\item{bootreps}{numeric-valued number of bootstrap replication for estimation
of 95\% confidence intervals. Default: 4000.}

\item{print_auc}{logical indicator for printing the area
under the ROC curve (AUC) on the plot.  Default: \verb{TRUE}.}

\item{partial_auc}{One of \verb{FALSE} or a length two numeric vector
of the form \verb{c(a, b)} where \verb{a} and \verb{b} are the endpoints
of the interval over which to compute the out-of-sample partial AUC (pAUC).
Ignored if \code{print_auc = FALSE}. Default: \code{c(0.8, 1)}.}

\item{partial_auc_focus}{one of \verb{"sensitivity"} or \verb{"specificity"},
indicating the measure for which the out-of-sample partial AUC is to be
computed. Default: \verb{"specificity"}.}

\item{partial_auc_correct}{logical indictor for transformation of the pAUC
to fall within the range from 0.5 (random guess) to 1.0 (perfect
classification). Default: \verb{TRUE}.}

\item{type}{type of plot; one of \verb{"l"} (line) or \verb{"S"}
(stair-step). Default: \verb{"l"}.}
}
\value{
This function produces a plot as a side effect, and (optionally)
returns a dataframe dataframe containing medians and
bootstrap confidence limits of sensitivity and specificity.
}
\description{
\code{plot.simple_screenr} plots the \emph{k}-fold cross-validated
receiver-operating characteristic, including confidence intervals on the
combinations of the local maxima of sensitivity and specificity.

Plot ROC curve with pointwise 95% confidence
intevals on sensitivity and specificity and (optionally) returns a dataframe
containing numerical values.
}
\examples{
data(unicorns)
too_simple <- simple_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
                          data = unicorns)
plot(too_simple)
}
\references{
Fawcett T. An introduction to ROC analysis. Pattern Recognition Letters. 2006.
27(8):861-874.
\url{https://doi.org/10.1016/j.patrec.2005.10.010}

Linden A. Measuring diagnostic and predictive accuracy in disease
management: an introduction to receiver operating characteristic (ROC) analysis.
Journal of Evaluation in Clinical Practice. 2006; 12(2):132-139.
\url{https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1365-2753.2005.00598.x}

Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez J-C, Muller M.
pROC: an open-source package for R and S+ to analyze and compare ROC curves.
BMC Bioinformatics 2011; 12:77. \url{https://www.biomedcentral.com/1471-2105/12/77}
}
