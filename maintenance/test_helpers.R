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
library(devtools)
library(dplyr)
#################################################################################
## Set paths and working directory
basepath <- file.path(Sys.getenv("DEVEL"), "screenr")
codepath <- file.path(basepath, "R")
workpath <- file.path(basepath, "maintenance")
datapath <- file.path(basepath, "data")
setwd(workpath)

#################################################################################
## Simulate loading the screenr package
#################################################################################
devtools::load_all()
attach(unicorns)

#################################################################################
## inverse_link
#################################################################################
new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
                        testresult = c(NA, NA), Q1 = c(0, 0), Q2 = c(0, 0),
                        Q3 = c(0, 0), Q4 = c(0, 0), Q5 = c(0, 1), Q6 = c(0, 1 ),
                        Q7 = c(0, 1))
mfit <- get_what(from = uniobj2, what = "ModelFit")
(coefs <- mfit$coefficients)
lp <- as.matrix(cbind(rep(1, nrow(new_corns)), new_corns[, 3:9])) %*%
           as.matrix(coefs, ncol =  1)
(preds <- inverse_link(lp, link = "logit"))
#################################################################################
## rescale_to_int
#################################################################################
x <- c(0.55, 1.21, 0.94, 0, 0.13)
rescale_to_int(x, max = 5)
#################################################################################
## roc_ci     NOTE: roc_ci is called from get_what( , what = "ROCci")
#################################################################################
(rocci <- get_what(from = uniobj2, what = "ROCci",  model = "minAIC" ))
debugonce(get_what.lasso_screenr)
debugonce(roc_ci )
aroc <- get_what(from = uniobj2, what = "cvROC",  model = "minAIC" )
(rocci <- roc_ci(aroc))
#################################################################################
## se_sp_max
#################################################################################
aroc <- aroc[, c("sensitivities", "specificities")]
aroc2 <- aroc
aroc2$specificities <- aroc2$specificities + 0.0001
aroc <- rbind(aroc, aroc2)
aroc <- aroc[order(aroc$sensitivities, decreasing = TRUE),]
se_sp_max(aroc)
#################################################################################
## sens_spec_plus
#################################################################################
(wtf <- with(unicorns, sens_spec_plus(test = Q1, gold = testresult,
                                      method = "jeffreys")))
class(wtf)
## $table
##              True +       True -      Total
## Test +           76         1371       1447
## Test -           21         4532       4553
## Total            97         5903       6000

## $ests
##                           est     lower     upper
## Apparent_positivity 0.2411667 0.2304688 0.2521141
## True_positivity     0.0161667 0.0132043 0.0195955
## Sensitivity         0.7835052 0.6939848 0.8564011
## Specificity         0.7677452 0.7568436 0.7783844
## PPV                 0.0525225 0.0419109 0.0649212
## NPV                 0.9953877 0.9930909 0.9970562

## R> 76 / (76 + 21)
## [1] 0.783505
## R> 4532 / (1371 + 4532)
## [1] 0.767745
## R> 76 / (76 + 1371)
## [1] 0.0525225
## R> 4532 / (4532 + 21)
## [1] 0.995388
debugonce(sens_spec_plus )
with(unicorns, sens_spec_plus(test = Q1, gold = testresult, method = "jeffreys"))
