

library(screenr)

?lasso_screenr
?logreg_screenr
?simple_screenr
?easy_tool


data(unicorns)
help(unicorns)
uniobj1 <- lasso_screenr(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7,
                          data = unicorns, Nfolds = 10)
summary(uniobj1)
print(uniobj1)

help(ntpp )
