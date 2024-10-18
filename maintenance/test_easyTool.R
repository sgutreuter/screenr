#################################################################################
##      R PROGRAM: test_easytool.R
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
library(pROC )
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


## uniobj1 <- lasso_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
##                          data = unicorns, Nfolds = 10, seed = 123)
## Test four cases:

wtf1 <- easy_tool(uniobj1, model = "minBIC", max = 3)
str(wtf1)
plot(wtf1)
summary(wtf1)
wtf2 <- easy_tool(uniobj1,  model =  "minBIC", max = 3, crossval =  FALSE)
wtf3 <- easy_tool(uniobj2, max = 3)
wtf4 <- easy_tool(uniobj2, max = 3, crossval = FALSE )

## Methods:
(get_what(from = wtf1, what = "QuestionWeights"))
(get_what(from = wtf2, what = "QuestionWeights"))
(get_what(from = wtf3, what = "QuestionWeights"))
(get_what(from = wtf4, what = "QuestionWeights"))
(get_what(from = wtf1, what = "Call"))

ntpp(wtf2)
ntpp(wtf3)

plot(wtf2)
plot(wtf3)

print(wtf2)

summary(wtf2)
summary(wtf3)

aroc <- get_what(from = wtf1, what = "ROC" )
roctbl <- get_what(from = wtf1, what = "ROCci", se.min = 0.5)
