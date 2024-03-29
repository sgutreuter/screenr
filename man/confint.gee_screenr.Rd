% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confint-methods.R
\name{confint.gee_screenr}
\alias{confint.gee_screenr}
\title{An S3 Method to Compute Confidence Limits from \code{gee_screenr} Objects}
\usage{
\method{confint}{gee_screenr}(
  object,
  ...,
  intercept = TRUE,
  or = FALSE,
  conf_level = 0.95,
  digits = 4
)
}
\arguments{
\item{object}{an object of class \code{gee_screenr}.}

\item{...}{optional arguments passed to \code{predict} methods.}

\item{intercept}{(logical) retain (\code{TRUE}, default) or drop
(\code{FALSE}) the intercept coefficients.}

\item{or}{return odds ratios if \verb{TRUE}. Default: FALSE (returns
logit-scale coefficients).}

\item{conf_level}{confidence level for normal-theory Wald-type confidence
intervals. Default: 0.95.}

\item{digits}{number of decimal places to be printed. Default: 4.}
}
\value{
\code{confint.gee_screenr} returns a dataframe containing the
estimated coefficients (or odds ratios) and their Wald-type
lower and upper confidence limits (lcl and ucl, respectively).
}
\description{
\code{confint.logreg_screenr} returns the logistic model parameter estimates
and their and Wald-type confidence limits from \code{gee_screenr}-class
objects.
}
\examples{
attach(uniobj3)
confint(uniobj3, or = TRUE)

}
