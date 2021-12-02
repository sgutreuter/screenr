#################################################################################
##      R PROGRAM: test_glmpathScreening.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: Testing sandbox for glmpathScreenr
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
source(file.path(codepath, "easyTool.R"))
source(file.path(codepath, "glmpathScreening.R"))
source(file.path(codepath, "logisticScreening.R"))
source(file.path(codepath, "helperFunctions.R"))
##load(file.path(datapath, "unicorns.rda"))
uniobj1 <- readRDS(file.path(datapath, "uniobj1.Rdata"))
uniobj2 <- readRDS(file.path(datapath, "uniobj2.Rdata"))

## Test four cases:

wtf1 <- easyTool(uniobj1, model = "minBIC", max = 3)
str(wtf1)

wtf2 <- easyTool(uniobj1, max = 3, crossval = FALSE )
str(wtf2)

wtf3 <- easyTool(uniobj2, max = 3)
str(wtf3)

wtf4 <- easyTool(uniobj2, max = 3, crossval = FALSE )
str(wtf4)

rm(wtf1, wtf2, wtf3, wtf4)

## Methods:
(getWhat(from = wtf1, what = "QuestionWeights"))
(getWhat(from = wtf1, what = "Call"))

ntpp(wtf2)
ntpp(wtf3)

plot(wtf2)
plot(wtf3)

print(wtf2)

summary(wtf2)
summary(wtf3)

aroc <- getWhat(from = wtf1, what = "ROC" )
