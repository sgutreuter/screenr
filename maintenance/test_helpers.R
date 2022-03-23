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
source(file.path(codepath, "helper_functions.R"))
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


roc <- get_what(from = uniobj1, what = "ROCci", se.min = 0.9)
print(roc)

## Test sens_spec_plus
sens_spec_plus(test = "Q1", gold = "testresult", data = unicorns)
## $table
##              True +       True -      Total
## Test +           76         1371       1447
## Test -           21         4532       4553
## Total            97         5903       6000

## $ests
##                           est     lower     upper
## Apparent_positivity 0.2411667 0.2303868 0.2521987
## True_positivity     0.0161667 0.0131291 0.0196868
## Sensitivity         0.7835052 0.6883042 0.8607100
## Specificity         0.7677452 0.7567575 0.7784677
## PPV                 0.0525225 0.0416017 0.0653015
## NPV                 0.9953877 0.9929581 0.9971427

## R> 76 / (76 + 21)
## [1] 0.783505
## R> 4532 / (1371 + 4532)
## [1] 0.767745
## R> 76 / (76 + 1371)
## [1] 0.0525225
## R> 4532 / (4532 + 21)
## [1] 0.995388
sens_spec_plus(test = "Q1", gold = "testresult", data = unicorns,
               method = "wilson")
