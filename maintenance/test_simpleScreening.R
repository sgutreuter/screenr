#################################################################################
##      R PROGRAM: simpleScreening_testing.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: Testing sandbox for simpleScreening
##
##     WRITTEN BY: Steve Gutreuter
##                 E-mail:  sgutreuter@gmail.gov
#################################################################################

#################################################################################
## Set paths and working directory
#################################################################################
codepath <- file.path(Sys.getenv("DEVEL"), "screenr/R")
workpath <- file.path(Sys.getenv("DEVEL"), "screenr/maintenance")
setwd(workpath)

library(pROC)
library(tidyverse)

#################################################################################
## Source the screenr R code for experimentation and testing
#################################################################################
source(file.path(codepath, "helperFunctions.R"))
source(file.path(codepath, "simpleScreening.R"))
load(file.path(datapath, "unicorns.rda"))
#################################################################################
## or attach the package
#################################################################################
##library(screenr)

#################################################################################
## Create a simpleScreenr object
#################################################################################
simp <- simpleScreenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6, data = unicorns)
str(simp)

#################################################################################
## Methods testing
#################################################################################

## getWhat
debugonce(getWhat.simpleScreenr )
roc <- getWhat(from = simp, what = "isROC" )
plot(roc)

## ntpp
ntpp(simp)

## plot
plot(simp)
plot(simp, plot_ci = FALSE)
plot(simp, print = FALSE)

## print
print(simp)

## summary
summary(simp )

#################################  End of File  #################################
