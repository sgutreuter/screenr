################################################################################
##      R PROGRAM: pROC_examples.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: Testing sandbox for the pROC package
##
##     WRITTEN BY: Steve Gutreuter
##                 E-mail:  sgutreuter@gmail.gov
################################################################################

## Examples from pROC

library(pROC)

data(aSAH)
## Not run:
## Start a ROC plot
rocobj <- roc(aSAH$outcome, aSAH$s100b)



plot(rocobj)
## Thresholds
ci.thresolds.obj <- ci.thresholds(rocobj, thresholds = "local maximas")
plot(ci.thresolds.obj)
## Sensitivities
plot(rocobj) # restart a new plot
ci.se.obj <- ci(rocobj, of="se", boot.n=2000)


## NOTE: The object returned by ci.se has everything needed to plot the ROC
##       with bootstrap CIs on sensitivity.
plot(ci.se.obj)
## Plotting a shape. We need more

ci.se.obj <- ci.se(rocobj, specificities=seq(0, 1, .01), boot.n=100)
plot(rocobj) # restart a new plot
plot(ci.se.obj, type="shape")
plot(ci.se.obj, type="bars")
## Direct syntax (response, predictor):
plot.roc(aSAH$outcome, aSAH$s100b,
ci=TRUE, of="thresholds")
