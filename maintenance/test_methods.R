#################################################################################
##      R PROGRAM: test_methods.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: S3 methods testing
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
attach(uniobj1.rda)
attach(uniobj2)
attach(val_data.rda)

#################################################################################
## coef methods
#################################################################################
(wtf1 <- coef(uniobj1))
(wtf2 <- coef(uniobj2))

lmobj <- get_what(uniobj2, what = "ModelFit")

coef(uniobj2)
coef(uniobj2, or = TRUE)
coef(uniobj2, or = TRUE, conf_level = 0.9)

#################################################################################
## plot methods
#################################################################################
## plot.easy_tool
debugonce(easy_tool)
tool <- easy_tool(uniobj1, model =  "minAIC", max = 3, crossval = TRUE)
plot(tool, type = "l")
plot(tool, print.auc = FALSE, plot_ci = FALSE, partial_auc =  FALSE)

## plot.lasso_creenr
plot(uniobj1, model = "minAIC", plot_ci =  FALSE, print_auc =  FALSE)

## plot.logreg_screenr
plot(uniobj2, print_auc = FALSE, plot_ci = FALSE)

## plot.simple_screenr
too_simple <- simple_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
                             data = unicorns)
plot(too_simple)

## plot.roc
new_preds <- predict(uniobj1, newdata = val_data)
new_roc <- pROC::roc(testresult ~ phat_minAIC, data = new_preds, auc = TRUE)
plot(new_roc,  print.auc = TRUE, partial.auc = c(0.8, 1.0), type = "S",
     partial.auc.focus = "sensitivity", partial.auc.correct = TRUE)

#################################################################################
## get_what methods
#################################################################################
## get_what.easyTool
tool <- easy_tool(uniobj2, max = 3, crossval = TRUE)
(get_what(from = tool,  what = "ROC") )
(get_what(from = tool,  what = "Call") )
(get_what(from = tool,  what = "Scores") )
(get_what(from = tool, what = "ROCci") )
(get_what(from = tool, what = "QuestionWeights") )
tool <- easy_tool(uniobj1, max = 3, crossval = TRUE)
(get_what(from = tool,  what = "ROC") )
(get_what(from = tool,  what = "Call") )
(get_what(from = tool,  what = "Scores") )
(get_what(from = tool, what = "ROCci") )
(get_what(from = tool, what = "QuestionWeights") )

## get_what.lasso_screenr methods
plot(get_what(from = uniobj1, what = "glmpathObj", model = "minAIC"))
plot(get_what(from = uniobj1, what = "glmpathObj", model = "minBIC"))
(get_what(from = uniobj1,  what = "ROCci", model = "minBIC"))

## getWhat.logreg_screenr
myROCci <- get_what(from = uniobj2, what = "ROCci")
print(myROCci)

## getWhat.simple_screenr
too_simple <- simple_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
                             data = unicorns)
myroc <- get_what(from = too_simple, what = "ROCci" )
print(myroc)

################################################################################
## Predict methods
################################################################################
newdata <- data.frame(ID = c("Alice D.", "Bernie P."),
                        testresult = c(NA, NA),
                        Q1 = c(0, 0), Q2 = c(0, 0), Q3 = c(0, 1), Q4 = c(0, 0),
                      Q5 = c(0, 1), Q6 = c(0, 1), Q7 = c(0, 1))
## predict.lasso_screenr
newpreds <- predict(uniobj1, newdata = newdata)

## predict.easy_tool
tool <- easy_tool(uniobj1, max = 3, crossval = TRUE)
print(predict(tool, newdata = newdata) )


#################################################################################
## ntpp methods
#################################################################################

se <- c(0.1, 0.5, 0.9)
sp <- se
prev <- rep(0.1,  3)
ntpp(se = se, sp = sp, prev = prev)
ntpp(tool)
(wtf <- ntpp(uniobj1))
ntpp(uniobj2)
obj0 <- data.frame(sensitivities = se,  specificities = sp,  prev = prev)
with(obj0, ntpp(se = sensitivities, sp =  specificities, prev =  prev))

## Test with se = 1 and sp = 0
ssp1m <- readRDS(file.path(workpath, "ssp1m.rds")) ## degenerate se and sp
(wtf <- with(ssp1m, ntpp(se = ests$est[3], sp = ests$est[4], prev =  ests$est[2] ) ))
