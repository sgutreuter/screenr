#################################################################################
##      R PROGRAM: test_glmpathScreening.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: Testing sandbox for binomialScreening
##
##     WRITTEN BY: Steve Gutreuter
##                 E-mail:  sgutreuter@gmail.gov
#################################################################################
library(stringr)
library(glmpath )
#################################################################################
## Set paths and working directory
codepath <- file.path(Sys.getenv("DEVEL"), "screenr/R")
workpath <- file.path(Sys.getenv("DEVEL"), "screenr/maintenance")
datapath <- file.path(Sys.getenv("DEVEL"), "screenr/data")
setwd(workpath)

#################################################################################
## Source the screenr R code and load the unicorns data
#################################################################################
source(file.path(codepath, "glmpathScreening.R"))
source(file.path(codepath, "helperFunctions.R"))
load(file.path(datapath, "unicorns.rda") )
uniobj <- readRDS(file.path(datapath, "uniobj.rds") )

#################################################################################
## Create new data for prediction
#################################################################################
new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
                        testresult = c(NA, NA), Q1 = c(0, 0),
                        Q2 = c(0, 0), Q3 = c(0, 1), Q4 = c(0, 0), Q5 = c(0, 1))


#################################################################################
## glmpathScreener testing
#################################################################################
uniobj <- glmpathScreener(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5,
                          data = unicorns, Nfolds = 2, seed = 123)
saveRDS(uniobj, file.path(datapath, "uniobj.rds" ))
print(uniobj)
summary(uniobj )
plot(uniobj, model =  "minAIC")
coef(uniobj, or = TRUE, intercept = FALSE)
coef(uniobj, or = FALSE, intercept = TRUE)

##debugonce(glmpath::predict.glmpath)
debugonce(predict.glmpathScreener )
(new_preds <- predict(uniobj, new_corns ))

#################################  End of File  #################################
