## Generate unicorn data
library(mvtnorm)
set.seed(321)
basepath <- file.path(Sys.getenv("DEVEL"), "screenr/R")
datapath <- file.path(Sys.getenv("DEVEL"), "screenr/data")
unimaker <- function(N, clusters = NULL, revar = NULL){
    nc <- length(clusters)
    if(!(N %% nc) == 0) stop("N must be an integer multiple of the number of clusters")
    reps <- as.integer(N / nc)
    q1 <- rbinom(N, 1, 0.1)
    q2 <- rbinom(N, 1, 0.05)
    q3 <- rbinom(N, 1, 0.15)
    q4 <- rbinom(N, 1, 0.05)
    q5 <- rbinom(N, 1, 0.02)
    parms <- c(-4.5, 2.0, 1.5, 1.5, 1.5, 2.8 )
    lp <- as.vector(matrix(c(rep(1, N), q1, q2, q3, q4, q5), ncol = 6) %*%
        matrix(parms, nrow = 6))
    cov <- matrix(rep(-0.01, nc^2), nrow = nc)
    diag(cov) <- 1
    covmat <- revar*cov
    remat <- as.vector(rmvnorm(1, mean = rep(0, nc), sigma = covmat))
    re <- NULL
    clus <- character()
    for(i in 1:nc){
        re <- c(re, rep(remat[i], reps))
        clus <- c(clus, rep(as.character(clusters[i]), reps))
    }
    p <- exp(lp + re) / ( 1 + exp(lp + re))
    y <- rbinom(N, 1, p)
    result <- data.frame(clinic = clus, Q1 = q1, Q2 = q2, Q3 = q3, Q4 = q4, Q5 = q5,
                         testresult = y)
}

clusters <- factor(1:20, labels = paste(rep("C", 20), as.character(1:20), sep="-"),
                   ordered = TRUE)
##revar = 0.01        ## original
revar = 0.05

unicorns <- unimaker(5000, clusters, revar)
rm(clusters, revar, unimaker)
##mean(unicorns$testresult)
##save(unicorns, file = file.path(workpath, "unicorns.rda"))
