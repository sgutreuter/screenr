#################################################################################
##      R PROGRAM: test_helpers.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: Testing sandbox for helper functions
##
##     WRITTEN BY: Steve Gutreuter
##                 E-mail:  sgutreuter@gmail.gov
#################################################################################
library(stringr)
library(glmpath)
#################################################################################
## Set paths and working directory
codepath <- file.path(Sys.getenv("DEVEL"), "screenr/R")
workpath <- file.path(Sys.getenv("DEVEL"), "screenr/maintenance")
datapath <- file.path(Sys.getenv("DEVEL"), "screenr/data")
setwd(workpath)

#################################################################################
## Source the screenr R code and load the unicorns data
#################################################################################
source(file.path(codepath, "helperFunctions.R"))
load(file.path(datapath, "unicorns.rda"))
load(file.path(datapath, "uniobj2.rda"))

#################################################################################
## Check inverseLink
#################################################################################

new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
                        testresult = c(NA, NA), Q1 = c(0, 0), Q2 = c(0, 0),
                        Q3 = c(0, 0), Q4 = c(0, 0), Q5 = c(0, 1), Q6 = c(0, 1 ),
                        Q7 = c(0, 1))
mfit <- getWhat(from = uniobj2, what = "ModelFit")
coefs <- mfit$coefficients
lp <- as.matrix(cbind(rep(1, nrow(new_corns)), new_corns[, 3:9])) %*%
           as.matrix(coefs, ncol =  1)
(preds <- inverseLink(lp, link = "logit"))
