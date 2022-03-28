## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
  )
library(screenr)
data(unicorns)
data(val_data)

## -----------------------------------------------------------------------------
## The first six lines of the unicorn screening data:
head(unicorns)

## -----------------------------------------------------------------------------
uniobj1 <- lasso_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
                         data = unicorns, Nfolds = 10, seed = 123)
class(uniobj1)

## -----------------------------------------------------------------------------
coef(uniobj1)

## -----------------------------------------------------------------------------
coef(uniobj1, or = TRUE, intercept = FALSE)

## ---- fig.width = 4, fig.height = 4, fig.fullwidth = TRUE---------------------
pathobj <- get_what(from = uniobj1, what = "glmpathObj", model = "minAIC")
plot(pathobj)

## -----------------------------------------------------------------------------
methods(class = "lasso_screenr")

## ---- eval = FALSE------------------------------------------------------------
#  saveRDS(uniobj1, file = "uniobj1.rds" )

## ---- fig.width = 4, fig.height = 4, fig.fullwidth = TRUE---------------------
plot(uniobj1, model = "minAIC")

## -----------------------------------------------------------------------------
roc_maximas <- get_what(from = uniobj1, what = "ROCci", se.min = 0.9)
print(roc_maximas)

## -----------------------------------------------------------------------------
new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
                        testresult = c(NA, NA),
                        Q1 = c(0, 0), Q2 = c(0, 0), Q3 = c(0, 1),
						Q4 = c(0, 0), Q5 = c(0, 1), Q6 = c(0, 0),
						Q7 = c(0, 0))
new <- predict(uniobj1,  newdata = new_corns )
print(new)

## -----------------------------------------------------------------------------
new_preds <- predict(uniobj1, newdata = val_data)
head(new_preds)

## -----------------------------------------------------------------------------
new_roc <- pROC::roc(testresult ~ phat_minAIC, data = new_preds,
auc = TRUE)
class(new_roc)

## ---- fig.width = 4, fig.height = 4-------------------------------------------
plot(new_roc, print.auc = TRUE, partial.auc =  c(0.8, 1.0),
     partial.auc.focus = "sensitivity", partial.auc.correct = TRUE)

## -----------------------------------------------------------------------------
new_perf <- roc_ci(new_roc, se.min = 0.8)
print(new_perf)

## -----------------------------------------------------------------------------
et_3 <- easy_tool(uniobj1, max = 3, model = "minAIC", crossval = TRUE)
class(et_3)

## -----------------------------------------------------------------------------
methods(class = "easy_tool")

## -----------------------------------------------------------------------------
qwts <- get_what(from = et_3, what = "QuestionWeights")
print(qwts)

## ---- fig.width = 4, fig.height = 4, fig.fullwidth = TRUE---------------------
plot(et_3)

## -----------------------------------------------------------------------------
qw_maximas <- get_what(from = et_3, what = "ROCci")
print(qw_maximas)

## ---- fig.width = 4, fig.height = 3-------------------------------------------
knitr::include_graphics("UniTool.png")

## -----------------------------------------------------------------------------
ntpp(et_3)

