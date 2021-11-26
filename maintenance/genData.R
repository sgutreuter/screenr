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
library(mvtnorm)
## Set paths and working directory
codepath <- file.path(Sys.getenv("DEVEL"), "screenr/R")
workpath <- file.path(Sys.getenv("DEVEL"), "screenr/maintenance")
datapath <- file.path(Sys.getenv("DEVEL"), "screenr/data")
setwd(workpath)
#################################################################################
## Function unimaker generates fake data
#################################################################################
unimaker <- function(N, seed = Sys.time( )){
    ## Generate the covariates
    q1 <- rbinom(N, 1, 0.1)
    q2 <- rbinom(N, 1, 0.05)
    q3 <- rbinom(N, 1, 0.02)
    q4 <- rbinom(N, 1, 0.15)
    q5 <- rbinom(N, 1, 0.05)
    q6 <- rbinom(N, 1, 0.02)
    ## Specify coefficient values
    parms <- c(-4.5, 2.0, 1.5, 0.05, 1.0, 1.5, 2.8)
    ## Compute the response probability and generate the response variable
    lp <- as.vector(matrix(c(rep(1, N), q1, q2, q3, q4, q5, q6), ncol = 7) %*%
        matrix(parms, nrow = 7))
    p <- exp(lp) / (1 + exp(lp))
    y <- rbinom(N, 1, p)
    result <- data.frame(Q1 = q1, Q2 = q2, Q3 = q3, Q4 = q4, Q5 = q5, Q6 = q6,
                         testresult = y)
}
#################################################################################
## Generate the data
#################################################################################
unicorns <- unimaker(5000, seed = 321)
mean(unicorns$testresult)
#################################################################################
## Save the data into the package data directory
#################################################################################
save(unicorns, file = file.path(datapath, "unicorns.rda"))
