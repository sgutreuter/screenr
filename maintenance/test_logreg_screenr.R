#################################################################################
##      R PROGRAM: test_gee_screenr.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: Testing sandbox for logreg_screenr
##
##     WRITTEN BY: Steve Gutreuter
##                 E-mail:  sgutreuter@gmail.gov
#################################################################################
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

#################################################################################
## Create and save a logreg_screenr object
#################################################################################
uniobj2 <- logreg_screenr(testresult ~ Q1 + Q2 + Q3 + Q5 + Q6 + Q7,
                          data = unicorns, link = "logit", Nfolds = 10)
save(uniobj2, file = file.path(datapath, "uniobj2.rda"), compress = "xz")
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

## confint
confint(uniobj2)
confint(uniobj2, or = TRUE)

## getWhat
mfit <- get_what(from = uniobj2, what = "ModelFit")
mfit$coef

## ntpp
ntpp(uniobj2)

## plot
plot(uniobj2)
plot(uniobj2, plot_ci = FALSE)
plot(uniobj2, print_auc = FALSE)

## predict
(preds <- predict(uniobj2, newdata = new_corns, type = "response"))

## print
print(uniobj2)

## summary
summary(uniobj2)

## easy_tool
(et2 <- easy_tool(uniobj2) )

#################################################################################
## Misc functions
#################################################################################
## inverseLink (helperFunctions.R)
mfit <- get_what(from = uniobj2, what = "ModelFit")
coefs <- mfit$coefficients
lp <- as.matrix(cbind(rep(1, nrow(new_corns)), new_corns[, 3:8])) %*%
            as.matrix(coefs, ncol =  1)
inverse_link(lp, link = "logit")


#################################  End of File  #################################
