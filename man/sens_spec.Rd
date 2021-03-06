% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helperFunctions.R
\name{sens_spec}
\alias{sens_spec}
\title{Compute Sensitivity and Specificity from a 2 x 2 Table}
\usage{
sens_spec(x)
}
\arguments{
\item{x}{a 2 x 2 table, with columns representing frequencies of
gold-standard status and rows representing frequencies of status ascertained
from testing.  The first row contains frequencies of negative test results
and the first column contain frequencies of true negatives.}
}
\value{
a list containing components sensitivity and specificity.
Sensitivities and specificities are displayed as proportions rather than
percentages.
}
\description{
Computes sensitivity and specificity of a test.
}
\examples{
Gold <- rbinom(20, 1, 0.50)
Test <- Gold; Test[c(3, 9, 12, 16)] <- 1 - Test[c(3, 9, 12, 16)]
sens_spec(table(Test, Gold))
}
