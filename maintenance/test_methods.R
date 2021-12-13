#################################################################################
##      R PROGRAM: test_methods.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: S3 methods testing
##
##     WRITTEN BY: Steve Gutreuter
##                 E-mail:  sgutreuter@gmail.gov
#################################################################################
library(stringr)
library(glmpath)
library(pROC )
#################################################################################
## Set paths and working directory
codepath <- file.path(Sys.getenv("DEVEL"), "screenr/R")
workpath <- file.path(Sys.getenv("DEVEL"), "screenr/maintenance")
datapath <- file.path(Sys.getenv("DEVEL"), "screenr/data")
setwd(workpath)

#################################################################################
## Source the screenr R code and load the unicorns data
#################################################################################
source(file.path(codepath, "lasso_screenr.R"))
source(file.path(codepath, "logreg_screenr.R"))
source(file.path(codepath, "simple_screenr.R"))
source(file.path(codepath, "easy_tool.R"))
source(file.path(codepath, "helper_functions.R"))
source(file.path(codepath, "coef-methods.R"))
source(file.path(codepath, "get_what-methods.R"))
source(file.path(codepath, "ntpp-methods.R"))
source(file.path(codepath, "predict-methods.R"))
source(file.path(codepath, "plot-methods.R"))
source(file.path(codepath, "print-methods.R"))
source(file.path(codepath, "summary-methods.R"))

load(file.path(datapath, "unicorns.rda"))
load(file.path(datapath, "uniobj1.rda"))
load(file.path(datapath, "uniobj2.rda"))
load(file.path(datapath, "val_data.rda"))

#################################################################################
## plot methods
#################################################################################
## plot.easy_tool
tool <- easy_tool(uniobj1, max = 3, crossval = TRUE)
plot(tool)

## plot.lasso_creenr
plot(uniobj1, model = "minAIC")

## plot.logreg_screenr
plot(uniobj2)

## plot.simple_screenr
too_simple <- simple_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
                             data = unicorns)
plot(too_simple)

#################################################################################
## get_what methods
#################################################################################
## get_what.easyTool
tool <- easy_tool(uniobj1, max = 3, crossval = TRUE)
ROCci <- get_what(from = tool, what = "ROCci")
print(ROCci)

## get_what.lasso_screenr methods
pathobj <- get_what(from = uniobj1, what = "glmpathObj", model = "minAIC")
plot(pathobj)
cvROCci <- get_what(from = uniobj1,  what = "ROCci", model = "minBIC")
print(cvROCci)

## getWhat.logreg_screenr
myROCci <- get_what(from = uniobj2, what = "ROCci")
print(myROCci)

## getWhat.simple_screenr
too_simple <- simple_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
                             data = unicorns)
myroc <- get_what(from = too_simple, what = "ROCci" )
print(myroc)


newpreds <- predict(uniobj1, newdata = val_data)
