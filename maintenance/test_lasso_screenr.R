#################################################################################
##      R PROGRAM: test_lasso_screenr.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: Testing sandbox for lasso_screenr
##
##     WRITTEN BY: Steve Gutreuter
##                 E-mail:  sgutreuter@gmail.gov
#################################################################################
library(stringr)
library(glmpath)
library(pROC)
#################################################################################
## Set paths and working directory
#################################################################################
codepath <- file.path(Sys.getenv("DEVEL"), "screenr/R")
workpath <- file.path(Sys.getenv("DEVEL"), "screenr/maintenance")
datapath <- file.path(Sys.getenv("DEVEL"), "screenr/data")
setwd(workpath)

#################################################################################
## Source the screenr R code and load the unicorns data
#################################################################################
source(file.path(codepath, "lasso_screenr.R"))
source(file.path(codepath, "coef-methods.R"))
source(file.path(codepath, "print-methods.R"))
source(file.path(codepath, "plot-methods.R"))
source(file.path(codepath, "summary-methods.R"))
source(file.path(codepath, "ntpp-methods.R"))
source(file.path(codepath, "predict-methods.R"))
source(file.path(codepath, "get_what-methods.R"))
source(file.path(codepath, "helper_functions.R"))
load(file.path(datapath, "unicorns.rda"))

#################################################################################
## Create and save a lasso_screenr object
#################################################################################
debugonce(lasso_screenr)
uniobj1 <- lasso_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
                         data = unicorns, Nfolds = 10, seed = 123)
##save(uniobj1, file = file.path(datapath, "uniobj1.rda"), compress = "xz")
summary(uniobj1)
coef(uniobj1)
coef(uniobj1, or = TRUE, intercept = FALSE)
plot(uniobj1)
pathobj <- get_what(from = uniobj1, what = "glmpathObj", model = "minAIC")
plot(pathobj)
print(uniobj1)

wtf <- get_what(from = uniobj1, what = "glmpathObj")
plot(wtf)
rm(wtf)

#################################################################################
## Create new data for prediction
#################################################################################
new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
                        testresult = c(NA, NA),
                        Q1 = c(0, 0), Q2 = c(0, 0), Q3 = c(0, 1), Q4 = c(0, 0),
                        Q5 = c(0, 1), Q6 = c(0, 0), Q7 = c(0, 0))

#################################################################################
## Methods testing
#################################################################################

## coef
coef(uniobj1)
coef(uniobj1, or = TRUE, intercept = FALSE)

## getWhat
pathobj <- get_what(from = uniobj1, what = "glmpathObj", model = "minAIC")
plot(pathobj)
cvROC <- get_what(from = uniobj1,  what = "cvROC", model = "minBIC")
plot(cvROC)

## ntpp
ntpp(uniobj1)

## plot
plot(uniobj1)
plot(uniobj1, print_ci = FALSE )

## predict
predict(uniobj1, newdata = new_corns)

## print
print(uniobj1)

## summary
summary(uniobj1)

debugonce(screenr:::predict.lasso_screenr)
new <- predict(uniobj1,  newdata = new_corns )
print(new )
#################################################################################
## Experimentation
#################################################################################

debugonce(plot.lasso_screenr)
plot(uniobj1)
plot(uniobj1, partial.auc =  FALSE )

#################################  End of File  #################################
