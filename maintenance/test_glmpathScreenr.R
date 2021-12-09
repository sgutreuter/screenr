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
source(file.path(codepath, "glmpathScreenr.R"))
source(file.path(codepath, "coef.R"))
source(file.path(codepath, "print.R"))
source(file.path(codepath, "summary.R"))
source(file.path(codepath, "getWhat.R"))
source(file.path(codepath, "helperFunctions.R"))
load(file.path(datapath, "unicorns.rda"))

#################################################################################
## Create and save a glmpathScreenr object
#################################################################################

uniobj1 <- glmpathScreenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
                          data = unicorns, Nfolds = 10, seed = 123)
##save(uniobj1, file = file.path(datapath, "uniobj1.rda"), compress = "xz")
coef(uniobj1)
print(uniobj1)
summary(uniobj1)
wtf <- getWhat(from = uniobj1, what = "glmpathObj")
plot(wtf)
rm(wtf)

#################################################################################
## Create new data for prediction
#################################################################################
new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
                        testresult = c(NA, NA),
                        Q1 = c(0, 0), Q2 = c(0, 0), Q3 = c(0, 1), Q4 = c(0, 0),
                        Q5 = c(0, 1), Q6 = c(0, 1), Q7 = c(0, 1))

#################################################################################
## Methods testing
#################################################################################

## coef
coef(uniobj1)
coef(uniobj1, or = TRUE, intercept = FALSE)

## getWhat
pathobj <- getWhat(from = uniobj1, what = "glmpathObj", model = "minAIC")
plot(pathobj)
cvROC <- getWhat(from = uniobj1,  what = "cvROC", model = "minBIC")
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

#################################################################################
## Experimentation
#################################################################################


#################################  End of File  #################################
