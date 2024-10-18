#################################################################################
##      R PROGRAM: genData.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: Generate the unicorn data
##
##     WRITTEN BY: Steve Gutreuter
##                 E-mail:  sgutreuter@gmail.gov
#################################################################################
library(tidyverse)
## Set paths and working directory
codepath <- file.path(Sys.getenv("DEVEL"), "screenr/R")
workpath <- file.path(Sys.getenv("DEVEL"), "screenr/maintenance")
datapath <- file.path(Sys.getenv("DEVEL"), "screenr/data")
setwd(workpath)
#################################################################################
## Function unimaker generates fake data
#################################################################################
unimaker <- function(N, parms = parms, seed = Sys.time( )){
    ## Generate the covariates
    set.seed(seed)
    q1 <- rbinom(N, 1, 0.24)
    q2 <- rbinom(N, 1, 0.04)
    q3 <- rbinom(N, 1, 0.06)
    q4 <- rbinom(N, 1, 0.23)
    q5 <- rbinom(N, 1, 0.20)
    q6 <- rbinom(N, 1, 0.10)
    q7 <- rbinom(N, 1, 0.71)
    ## Compute the response probability and generate the response variable
    lp <- as.vector(matrix(c(rep(1, N), q1, q2, q3, q4, q5, q6, q7), ncol = 8) %*%
        matrix(parms, nrow = 8))
    p <- exp(lp) / (1 + exp(lp))
    set.seed(seed + 1)
    y <- rbinom(N, 1, p)
    ID <- 1:length(y)
    result <- data.frame(ID =  ID, Q1 = q1, Q2 = q2, Q3 = q3, Q4 = q4, Q5 = q5,
                         Q6 = q6, Q7 = q7, testresult = y)
}
#################################################################################
## Generate and insert the training data
#################################################################################
parms1 <- c(-7.5, 2.5, 1.5, 2.75, 0.01, 1.5, 0.75, 1.2)
unicorns <- unimaker(N = 6000, parms =  parms1, seed = 999)
mean(unicorns$testresult)
save(unicorns, file = file.path(datapath, "unicorns.rda"))
#################################################################################
## Generate and insert validation data
#################################################################################
parms2 <- c(-7.3, 2.6, 1.4, 2.7, 0.01, 1.4, 0.85, 1.2 )
val_data <- unimaker(N = 3000, parms = parms2, seed = 333)
mean(val_data$testresult)
save(val_data, file = file.path(datapath, "val_data.rda"))
#################################################################################
