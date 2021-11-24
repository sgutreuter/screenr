#################################################################################
##      R PROGRAM: test_logisticScreening.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: Testing sandbox for binomialScreening
##
##     WRITTEN BY: Steve Gutreuter
##                 E-mail:  sgutreuter@gmail.gov
#################################################################################

#################################################################################
## Set paths and working directory
codepath <- file.path(Sys.getenv("DEVEL"), "screenr/R")
workpath <- file.path(Sys.getenv("DEVEL"), "screenr-testing")
datapath <- file.path(Sys.getenv("DEVEL"), "screenr/data")
setwd(workpath)

#################################################################################
## Source the screenr R code for experimentation and testing
#################################################################################
source(file.path(codepath, "logisticScreening.R"))
source(file.path(codepath, "helperFunctions.R"))
load(file.path(datapath, "unicorns.rda") )
load(file.path(datapath, "uniobj.Rdata") )

#################################################################################
## Create new data for prediction
#################################################################################
new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
                        testresult = c(NA, NA), Q1 = c(0, 0),
                        Q2 = c(0, 0), Q3 = c(0, 1), Q4 = c(0, 0), Q5 = c(0, 1))




bsobj2 <- logisticScreenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5,
                           data = unicorns, link = "logit", Nfolds = 10)

str(bsobj2)

summary(bsobj2)

print(bsobj2)

getROC(bsobj)
getROC(bsobj, simplify = FALSE)

debugonce(plot.logisticscreenr)
plot(bsobj2)
plot(bsobj2, print_ci = FALSE)

testCounts(bsobj2)

## Compute point estimates of their predicted probabilities testing positive:
predict(bsobj2$ModelFit, newdata = new, type = "response")
## or, compute directly
inverseLink("logit",
            as.matrix(cbind(rep(1, nrow(new)), new[, 2:6])) %*%
                            as.matrix(bsobj$ParamEst, ncol = 1))


#################################  End of File  #################################
