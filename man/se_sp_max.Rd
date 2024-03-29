% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{se_sp_max}
\alias{se_sp_max}
\title{Return a Simplified Dataframe of Sensitivity and Specificity}
\usage{
se_sp_max(object)
}
\arguments{
\item{object}{a dataframe containing at least columns named \code{sensitivities}
and \code{specificities}}
}
\value{
\code{se_sp_max} returns a dataframe which is a subset of \code{object}
containing only those rows for which specificity was the maximum for each
unique value of sensitivity.
}
\description{
Given a dataframe containing multiple values of specificity for
each value of sensitivity, return only the rows containing the largest value of
specificity for each unique value of sensitivity.
}
