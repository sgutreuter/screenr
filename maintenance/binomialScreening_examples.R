data(unicorns)
help(unicorns)
unitool <- binomialScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5,
                             data = unicorns, link = "logit", Nfolds = 20)
summary(unitool)
plot(unitool)
testCounts(unitool)

## Example implementation of screening based on those results
## Suppose there are new observations (excluding testing) from two previously
## untested unicorns:

new <- data.frame(ID = c('"Bernie P."', '"Alice D."'), Q1 = c(0, 0), Q2 = c(0, 0),
                    Q3 = c(1, 0), Q4 = c(0, 0), Q5 = c(1, 0))
print(new)

## Compute point estimates of their predicted probabilities testing positive:
inverseLink(as.matrix(cbind(rep(1, nrow(new)), new[, 2:6])) %*%
                           as.matrix(unitool$ParamEst, ncol = 1), "logit")
## or, more simply,
predict(unitool$ModelFit, newdata = new, type = "response")
