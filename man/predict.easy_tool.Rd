% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict-methods.R
\name{predict.easy_tool}
\alias{predict.easy_tool}
\title{An S3 Method to Compute Simplified Screening Scores}
\usage{
\method{predict}{easy_tool}(object = NULL, ..., newdata = NULL)
}
\arguments{
\item{object}{an object of class \code{easy_tool} produced by
\code{`easy_tool`}.}

\item{...}{optional arguments to \code{predict} methods.}

\item{newdata}{new dataframe from which predicted simplified screening
scores are desired. The dataframe must contain values of the same
response variables and covariates that were used to obtain \code{object}.}
}
\value{
\code{predict.easy_tool} returns (invisibly) a dataframe
augmenting \code{newdata} with the predicted simplified test screening
scores \code{score}.
}
\description{
\code{predict.easy_tool} computes predicted simplified
screening scores from new data.
}
\examples{
new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
                        testresult = c(NA, NA), Q1 = c(0, 0), Q2 = c(0, 0),
                        Q3 = c(0, 0), Q4 = c(0, 0), Q5 = c(0, 1), Q6 = c(0, 1 ),
                        Q7 = c(0, 1))
attach(uniobj1)
et <- easy_tool(uniobj1, max = 3)
print(predict(et, newdata = new_corns))
}
