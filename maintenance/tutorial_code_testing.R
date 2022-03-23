#################################################################################
##      R PROGRAM: tutorial_code_testing.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: Testing sandbox for lasso_screenr
##
##     WRITTEN BY: Steve Gutreuter
##                 E-mail:  sgutreuter@gmail.gov
#################################################################################
library(screenr)

workpath <- file.path(Sys.getenv("DEVEL"), "screenr/maintenance")
setwd(workpath)

data(unicorns)

uniobj1 <- lasso_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
                         data = unicorns, Nfolds = 10, seed = 123)

coef(uniobj1)
coef(uniobj1, or = TRUE, intercept = FALSE)

pathobj <- get_what(from = uniobj1, what = "glmpathObj", model = "minAIC")
plot(pathobj)

methods(class = "lasso_screenr")

plot(uniobj1, model = "minAIC")
roc_maximas <- get_what(from = uniobj1, what = "ROCci", se.min = 0.9)
print(roc_maximas)

new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
                        testresult = c(NA, NA),
                        Q1 = c(0, 0), Q2 = c(0, 0), Q3 = c(0, 1),
                        Q4 = c(0, 0), Q5 = c(0, 1), Q6 = c(0, 0),
                        Q7 = c(0, 0))
new <- predict(uniobj1,  newdata = new_corns )
print(new)

et_3 <- easy_tool(uniobj1, max = 3, model = "minAIC", crossval = TRUE)
class(et_3)
methods(class = "easy_tool")

qwts <- get_what(from = et_3, what = "QuestionWeights")
print(qwts)

plot(et_3)

qw_maximas <- get_what(from = et_3, what = "ROCci")
print(qw_maximas)

ntpp(et_3)

new_preds <- predict(uniobj1, newdata = val_data)
head(new_preds)
new_roc <- pROC::roc(testresult ~ phat_minAIC, data = new_preds, auc = TRUE)
plot(new_roc, print.auc = TRUE)

debugonce(roc_ci)
new_perf <- roc_ci(new_roc, se.min = 0.8, bootreps =  2000, conf.level = 0.95)
print(new_perf, conf.level =  0.95)
