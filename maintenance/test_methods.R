#################################################################################
##      R PROGRAM: test_methods.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: Testing sandbox for methods
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
source(file.path(codepath, "easyTool.R"))
source(file.path(codepath, "glmpathScreenr.R"))
source(file.path(codepath, "plot.R"))
source(file.path(codepath, "coef.R"))
source(file.path(codepath, "print.R"))
source(file.path(codepath, "predict.R"))
source(file.path(codepath, "summary.R"))
source(file.path(codepath, "getWhat.R"))
source(file.path(codepath, "ntpp.R"))
source(file.path(codepath, "simpleScreenr.R"))
source(file.path(codepath, "helperFunctions.R"))
load(file.path(datapath, "unicorns.rda"))
load(file.path(datapath, "uniobj1.rda"))
load(file.path(datapath, "uniobj2.rda"))

#################################################################################
## plot methods
#################################################################################
## plot.easyTool
tool <- easyTool(uniobj1, max = 3, crossval = TRUE)
plot(tool)

## plot.glmpathScreenr
plot(uniobj1, model = "minAIC")

## plot.logisticScreenr
plot(uniobj2)

## plot.simpleScreenr
too_simple <- simpleScreenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
                          data = unicorns)
plot(too_simple)

#################################################################################
## getWhat methods
#################################################################################
## getWhat.easyTool
tool <- easyTool(uniobj1, max = 3, crossval = TRUE)
ROCci <- getWhat(from = tool, what = "ROCci")
print(ROCci)
## TODO: Fix threshold variable name and round up values

## getWhat.glmpathScreenr methods
pathobj <- getWhat(from = uniobj1, what = "glmpathObj", model = "minAIC")
plot(pathobj)
cvROCci <- getWhat(from = uniobj1,  what = "ROCci", model = "minBIC")
print(cvROCci)

## getWhat.logisticScreenr
myROCci <- getWhat(from = uniobj2, what = "ROCci")
print(myROCci)

## getWhat.simpleScreenr
too_simple <- simpleScreenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6,
                            data = unicorns)
myroc <- getWhat(from = too_simple, what = "ROCci" )
print(myroc)
