#################################################################################
##      R PROGRAM: test_binomialScreening.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: Testing sandbox for binomialScreening
##
##     WRITTEN BY: Steve Gutreuter
##                 E-mail:  sgutreuter@gmail.gov
#################################################################################

#################################################################################
## Set paths and working directory
codepath <- file.path(Sys.getenv("DEVEL"), "screenr/R")
workpath <- file.path(Sys.getenv("DEVEL"), "screenr-testing")
setwd(workpath)

#################################################################################
## Source the screenr R code for experimentation and testing
#################################################################################
source(file.path(codepath, "binomialScreening.R"))
source(file.path(codepath, "helperFunctions.R"))
source(file.path(workpath, "genData.R"))

#################################################################################
## Create new data for prediction
#################################################################################
new <- data.frame(ID = c('"Bernie P."', '"Alice D."'), Q1 = c(0, 0), Q2 = c(0, 0),
                  Q3 = c(1, 0), Q4 = c(0, 0), Q5 = c(1, 0))
##print(new)

#################################################################################
## binomScreening
#################################################################################
binomialScreening <- function(formula,
                              data = NULL,
                              link = "logit",
                              Nfolds = 40L,
                              ...){
    if(!inherits(formula, "formula")) stop("Specify an model formula")
    if(!is.data.frame(data)) stop("Provide a data frame")
    if(!link %in% c("logit", "cloglog", "probit")) stop("Invalid link")
    call <- match.call()
    m <- match(c("formula", "data"), names(call), 0L)
    mf <- call[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    dat <- eval(mf, parent.frame())
    dat <- dat[complete.cases(dat), ]
    if(Nfolds > 0.20*dim(dat)[1])
        stop("Nfolds must be < 20% of number of complete observations")
    y <- stats::model.response(dat, "numeric")
    if(!all(y %in% c(0, 1))) stop("Response variable must be binary (0, 1)")
    prev <- mean(y, na.rm = TRUE)
    lrfit <- stats::glm(formula, data = dat, family = binomial(link = link))
    is.roc <- pROC::roc(lrfit$y, lrfit$fitted.values)
    N <- nrow(dat)
    holdouts <- split(sample(1:N), 1:Nfolds)
    cv.results <- data.frame(NULL)
    for(i in 1:Nfolds){
        res <- stats::glm(formula, data = dat[-holdouts[[i]], ],
                          family = binomial(link = link))
        pred.prob <- inverseLink(link, stats::predict(res, newdata = dat[holdouts[[i]], ]))
        y <- stats::model.response(dat[holdouts[[i]], ])
        cv.results <- rbind(cv.results,
                            data.frame(cbind(fold = rep(i, length(pred.prob)),
                                             y = y,
                                             data = dat[holdouts[[i]],],
                                             cv.pred.prob = pred.prob)))
    }
    cv.roc <- pROC::roc(cv.results$y, cv.results$cv.pred.prob,
                                  auc = TRUE)
    class(cv.results) <-  c("cv.predictions", "data.frame")
    result <- list(Call = call,
                   ModelFit = lrfit,
                   Prevalence = prev,
                   ParamEst = lrfit$coeff,
                   ISroc = is.roc,
                   CVpreds = cv.results,
                   CVroc = cv.roc)
    class(result) <- "binomscreenr"
    invisible(result)
}

plot.binomscreenr <- function(x, plot_ci = TRUE, conf_level = 0.95, bootreps = 2000,
                              print_ci = TRUE, ...){
    if(!class(x) == "binomscreenr") stop("x is not a binomscreenr object")
    stopifnot(conf_level > 0 & conf_level < 1)
    plot(x$CVroc, print.auc = TRUE, ci = FALSE, ...)
    if(plot_ci | print_ci){
        ciplt <- ci.thresholds(x$CVroc, boot.n = bootreps, progress = "text",
                               conf.level = conf_level, thresholds = "local maximas")
    }
    if(print_ci){
        threshold <- attr(ciplt, "thresholds")
        citable <- data.frame(cbind(threshold, ciplt$sensitivity,
                                    ciplt$specificity))
        names(citable) <- c("threshold", "se.low", "se.median",
                            "se.high", "sp.low", "sp.median", "sp.high")
        row.names(citable) <- 1:(dim(ciplt$sensitivity)[1])
    }
    if(plot_ci) plot(ciplt)
    lines.roc(x$ISroc, lty = 3)
    legend("bottomright", legend = c("cross-validated", "in-sample"), lty = c(1, 3),
           lwd = c(2, 2))
    if(print_ci) return(citable)
}


bsobj <- binomialScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5,
                           data = unicorns, link = "logit", Nfolds = 20)

str(bsobj)

summary(bsobj)

print(bsobj)

getROC(bsobj)
getROC(bsobj, simplify = FALSE)

debugonce(plot.binomscreenr)
plot(bsobj)
plot(bsobj, print_ci = FALSE)

testCounts(bsobj)

## Compute point estimates of their predicted probabilities testing positive:
predict(bsobj$ModelFit, newdata = new, type = "response")
## or, compute directly
inverseLink("logit",
            as.matrix(cbind(rep(1, nrow(new)), new[, 2:6])) %*%
                            as.matrix(bsobj$ParamEst, ncol = 1))


#################################  End of File  #################################
