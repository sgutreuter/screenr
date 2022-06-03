## ---- setup, include = FALSE--------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
  )
library(screenr)
data(unicorns)
data(val_data)

## ---- showdata----------------------------------------------------------------
## The first six lines of the unicorn screening data:
head(unicorns)

## ---- lasso_screenr-----------------------------------------------------------
uniobj1 <- lasso_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
                         data = unicorns, Nfolds = 10, seed = 123)
class(uniobj1)

## ---- coef-1------------------------------------------------------------------
coef(uniobj1)

## ---- coef-2------------------------------------------------------------------
coef(uniobj1, or = TRUE, intercept = FALSE)

## ---- pathplot, fig.width = 4, fig.height = 4, fig.fullwidth = TRUE-----------
pathobj <- get_what(from = uniobj1, what = "glmpathObj", model = "minAIC")
plot(pathobj)

## ---- class-lasso-------------------------------------------------------------
methods(class = "lasso_screenr")

## ---- save-rds, eval = FALSE--------------------------------------------------
#  saveRDS(uniobj1, file = "uniobj1.rds")

## ---- logreg-fit--------------------------------------------------------------
uniobj2 <- logreg_screenr(testresult ~ Q1 + Q2 + Q3 + Q5 + Q6 + Q7,
                         data = unicorns, Nfolds = 10, seed =  123)

## ---- coef-uniobj2------------------------------------------------------------
confint(uniobj2)

## ---- or-uniobj2--------------------------------------------------------------
confint(uniobj2, or = TRUE, intercept = FALSE)

## ---- lasso-plot, fig.width = 4, fig.height = 4, fig.fullwidth = TRUE---------
plot(uniobj2)

## ---- roc-maximas-lasso-------------------------------------------------------
roc_maximas <- get_what(from = uniobj2, what = "ROCci", se_min = 0.9)
print(roc_maximas)

## ---- new-corns-df------------------------------------------------------------
new_corns <- data.frame(ID = c("Alice D.", "Bernie P."),
                        testresult = c(NA, NA),
                        Q1 = c(0, 0), Q2 = c(0, 0), Q3 = c(0, 1),
					    Q5 = c(0, 1), Q6 = c(0, 0), Q7 = c(0, 0))
new <- predict(uniobj2, newdata = new_corns)
print(new)

## ---- get-val_data, results = "hide"------------------------------------------
val_data <- val_data[, -5]

## ---- new-preds---------------------------------------------------------------
new_preds <- predict(uniobj2, newdata = val_data)
head(new_preds)

## ---- new-roc-----------------------------------------------------------------
new_roc <- pROC::roc(testresult ~ phat, data = new_preds, auc = TRUE)
class(new_roc)

## ---- plot-new-roc, fig.width = 4, fig.height = 4-----------------------------
plot(new_roc, print_auc = TRUE, partial_auc =  c(0.8, 1.0), type = "S",
     partial_auc_focus = "sensitivity", partial_auc_correct = TRUE)

## ---- new-perf----------------------------------------------------------------
new_perf <- roc_ci(new_roc, se_min = 0.8)
print(new_perf)

## ---- make-et3----------------------------------------------------------------
et_3 <- easy_tool(uniobj2, max = 3, model = "minAIC", crossval = TRUE)
class(et_3)

## ---- easy-tool-methods-------------------------------------------------------
methods(class = "easy_tool")

## -----------------------------------------------------------------------------
qwts <- get_what(from = et_3, what = "QuestionWeights")
print(qwts)

## ---- plot-et3, fig.width = 4, fig.height = 4, fig.fullwidth = TRUE-----------
plot(et_3)

## ---- get_et3_maximas---------------------------------------------------------
qw_maximas <- get_what(from = et_3, what = "ROC")
print(qw_maximas)

## ---- show-simple-example, fig.width = 4, fig.height = 3----------------------
knitr::include_graphics("UniTool.png")

## ---- ntpp-et3----------------------------------------------------------------
ntpp(et_3)

