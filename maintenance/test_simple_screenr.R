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
library(devtools)

################################################################################
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
## or attach the package
#################################################################################
##library(screenr)

#################################################################################
## Create a simpleScreenr object
#################################################################################
simp <- simple_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6, data = unicorns)
str(simp)

#################################################################################
## Methods testing
#################################################################################

## getWhat
roc <- get_what(from = simp, what = "isROC" )
plot(roc)

## ntpp
ntpp(simp)

## plot
plot(simp)
plot(simp, plot_ci = FALSE)
plot(simp, print_auc = FALSE, plot_ci = FALSE)

## print
print(simp)

## summary
summary(simp )

#################################  End of File  #################################
