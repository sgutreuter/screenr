#################################################################################
##      R PROGRAM: test_logisticScreening.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: Testing sandbox for logreg_screenr
##
##     WRITTEN BY: Steve Gutreuter
##                 E-mail:  sgutreuter@gmail.gov
#################################################################################

#################################################################################
## Set paths and working directory
pkgpath <- file.path(Sys.getenv("DEVEL"), "screenr" )
codepath <- file.path(pkgpath, "R")
workpath <- file.path(pkgpath, "maintenance")
datapath <- file.path(pkgpath, "data")
setwd(workpath)

#################################################################################
## Source the screenr R code for experimentation and testing
#################################################################################
source(file.path(codepath, "logreg_screenr.R"))
source(file.path(codepath, "coef-methods.R"))
source(file.path(codepath, "print-methods.R"))
source(file.path(codepath, "summary-methods.R"))
source(file.path(codepath, "ntpp-methods.R"))
source(file.path(codepath, "predict-methods.R"))
source(file.path(codepath, "get_what-methods.R"))
source(file.path(codepath, "helper_functions.R"))
load(file.path(datapath, "unicorns.rda"))
##uniobj2 <- readRDS(file.path(datapath, "uniobj2.Rdata") )

#################################################################################
## Create and save a logreg_screenr object
#################################################################################
uniobj2 <- logreg_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
                          data = unicorns, link = "logit", Nfolds = 10)
##save(uniobj2, file = file.path(datapath, "uniobj2.rda"), compress = "xz")
coef(uniobj2)

#################################################################################
## Create new data for prediction
################################################################################
new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
                        testresult = c(NA, NA), Q1 = c(0, 0), Q2 = c(0, 0),
                        Q3 = c(0, 0), Q4 = c(0, 0), Q5 = c(0, 1), Q6 = c(0, 1 ))

#################################################################################
## Methods testing
#################################################################################

## coef
coef(uniobj2)
coef(uniobj2, intercept = FALSE, or = TRUE)

## getWhat
mfit <- get_what(from = uniobj2, what = "ModelFit")
mfit$coef

## ntpp
ntpp(uniobj2)

## plot
plot(uniobj2)
plot(uniobj2, print_ci = FALSE)

## predict
(preds <- predict(uniobj2, newdata = new_corns, type = "response"))

## print
print(uniobj2)

## summary
summary(uniobj2)

#################################################################################
## Misc functions
#################################################################################
## inverseLink (helperFunctions.R)
mfit <- getWhat.logreg_screenr(from = uniobj2, what = "ModelFit")
mfit <- getWhat(from = uniobj2, what = "ModelFit")
coefs <- mfit$coefficients
lp <- as.matrix(cbind(rep(1, nrow(new_corns)), new_corns[, 3:8])) %*%
            as.matrix(coefs, ncol =  1)
inverseLink(lp, link = "logit")


#################################  End of File  #################################
