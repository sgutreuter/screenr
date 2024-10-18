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
library(devtools)
#################################################################################
## Set paths and working directory
#################################################################################
codepath <- file.path(Sys.getenv("DEVEL"), "screenr/R")
workpath <- file.path(Sys.getenv("DEVEL"), "screenr/maintenance")
datapath <- file.path(Sys.getenv("DEVEL"), "screenr/data")
setwd(workpath)

#################################################################################
## Simulate loading the screenr package
#################################################################################
devtools::load_all()
attach(unicorns)
#################################################################################
## Create new data for prediction
#################################################################################
new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
                        testresult = c(NA, NA),
                        Q1 = c(0, 0), Q2 = c(0, 0), Q3 = c(0, 1), Q4 = c(0, 0),
                        Q5 = c(0, 1), Q6 = c(0, 0), Q7 = c(0, 0))
#################################################################################
## Create and save a lasso_screenr object
#################################################################################
##debugonce(lasso_screenr)
uniobj1 <- lasso_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
                         data = unicorns, Nfolds = 10, seed = 123)
##save(uniobj1, file = file.path(datapath, "uniobj1.rda"), compress = "xz")
#################################################################################
## Methods testing
#################################################################################
## summary
summary(uniobj1)
## coef
coef(uniobj1)
coef(uniobj1, or = TRUE, intercept = FALSE)
## getWhat
pathobj <- get_what(from = uniobj1, what = "glmpathObj", model = "minAIC")
plot(pathobj, xvar = "lambda")
cvROC <- get_what(from = uniobj1, what = "cvROC", model = "minBIC")
plot(cvROC)
## ntpp
ntpp(uniobj1)
ntpp(uniobj1, model = "minBIC")
## plot
plot(uniobj1)
plot(uniobj1, plot_ci = FALSE)
plot(uniobj1, model = "minBIC")
plot(uniobj1, model = "minBIC", type = "S")
## get and print the cross-validated ROC
roc_maximas <- get_what(from = uniobj1, what = "ROCci", se.min = 0.9)
print(roc_maximas)
## predict
predict(uniobj1, newdata = new_corns)
## print
print(uniobj1)
#################################  End of File  #################################
