

library(glmnet)

## Set paths and working directory
codepath <- file.path(Sys.getenv("DEVEL"), "screenr/R")
workpath <- file.path(Sys.getenv("DEVEL"), "screenr/maintenance")
datapath <- file.path(Sys.getenv("DEVEL"), "screenr/data")
setwd(workpath)

load(file.path(datapath, "unicorns.Rdata"))

obj <- glmnet(as.matrix(unicorns[, 2:7]), unicorns$testresult,
              family = binomial, standardize = FALSE)
coef(obj)



## R> coef(uniobj1)
##           AIC-best model BIC-best model
## Intercept      -5.001381      -5.001381
## Q1              2.410133       2.410133
## Q2              1.068597       1.068597
## Q3              0.264438       0.264438
## Q4              1.292780       1.292780
## Q5              1.525406       1.525406
## Q6              2.994698       2.994698
