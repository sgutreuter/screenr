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
    q1 <- rbinom(N, 1, 0.1)
    q2 <- rbinom(N, 1, 0.05)
    q3 <- rbinom(N, 1, 0.07)  ## 0.02
    q4 <- rbinom(N, 1, 0.15)
    q5 <- rbinom(N, 1, 0.05)
    q6 <- rbinom(N, 1, 0.02)
    ## Compute the response probability and generate the response variable
    lp <- as.vector(matrix(c(rep(1, N), q1, q2, q3, q4, q5, q6), ncol = 7) %*%
        matrix(parms, nrow = 7))
    p <- exp(lp) / (1 + exp(lp))
    set.seed(seed)
    y <- rbinom(N, 1, p)
    ID <- 1:length(y)
    result <- data.frame(ID =  ID, Q1 = q1, Q2 = q2, Q3 = q3, Q4 = q4, Q5 = q5,
                         Q6 = q6, testresult = y)
}
#################################################################################
## Generate and insert the training data
#################################################################################
parms1 <- c(-4.5, 2.0, 1.5, 0.0, 1.0, 1.5, 2.8)
unicorns <- unimaker(N = 6000, parms =  parms1, seed = 321)
mean(unicorns$testresult)
save(unicorns, file = file.path(datapath, "unicorns.rda"))
#################################################################################
## Generate and insert validation data
#################################################################################
parms2 <- c(-4.5, 2.01, 1.499, 0.0, 1.001, 1.501, 2.798)
val_data <- unimaker(N = 3000, parms = parms2, seed = 137)
mean(val_data$testresult )
save(val_data, file = file.path(datapath, "val_data.rda"))
