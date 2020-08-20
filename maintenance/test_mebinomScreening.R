#################################################################################
##      R PROGRAM: test_mebinomScreening.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: Testing sandbox for mebinomScreening
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
source(file.path(codepath, "helperFunctions.R"))
source(file.path(codepath, "mebinomScreening.R"))
source(file.path(codepath, "binomialScreening.R"))
source(file.path(workpath, "genData.R"))

#################################################################################
## Create new data for prediction
#################################################################################
new <- data.frame(ID = c('"Bernie P."', '"Alice D."'), Q1 = c(0, 0), Q2 = c(0, 0),
                  Q3 = c(1, 0), Q4 = c(0, 0), Q5 = c(1, 0))

#################################################################################
## mebinomScreening
#################################################################################
mebinomScreening <- function(formula,
                             id = NULL,
                             data = NULL,
                             link = "logit",
                             Nfolds = 40L,
                              ...){
    if(!inherits(formula, "formula")) stop("Specify an model formula")
    if(!is.data.frame(data)) stop("Provide a data frame")
    if(!link %in% c("logit", "cloglog", "probit")) stop("Invalid link")
    call <- match.call()
    meform <- update(formula, paste("~ . + (1|", id, ")"))
    formx <- update(formula, paste("~ . +", id))
    mf <- stats::model.frame(formx, data)
    dat <- eval(mf, parent.frame())
    dat <- dat[complete.cases(dat), ]
    if(Nfolds > 0.20*dim(dat)[1])
        stop("Nfolds must be < 20% of number of complete observations")
    y <- stats::model.response(dat, "numeric")
    if(!all(y %in% c(0, 1))) stop("Response variable must be binary (0, 1)")
    prev <- mean(y, na.rm = TRUE)
    lrfit <- lme4::glmer(meform, data = dat, family = binomial(link = link))
    is.roc <- pROC::roc(lrfit@resp$y, fitted(lrfit))
    N <- nrow(dat)
    holdouts <- split(sample(1:N), 1:Nfolds)
    cv.results <- data.frame(NULL)
    for(i in 1:Nfolds){
        res <- lme4::glmer(meform, data = dat[-holdouts[[i]], ],
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
                   ParamEst = lrfit@beta,
                   ISroc = is.roc,
                   CVpreds = cv.results,
                   CVroc = cv.roc)
    class(result) <- "binomscreenr"
    invisible(result)
}


meobj <- mebinomScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5, id = "clinic",
                          data = unicorns, link = "logit")

str(meobj)

summary(meobj)

print(meobj)

getROC(meobj)
getROC(meobj, simplify = FALSE)

plot(meobj)
plot(meobj, print_ci = FALSE)

## Compute point estimates of their predicted probabilities testing positive:
predict(meobj$ModelFit, newdata = new, type = "response", re.form = NA)
## or, compute directly
inverseLink("logit",
            as.matrix(cbind(rep(1, nrow(new)), new[, 2:6])) %*%
                            as.matrix(meobj$ParamEst, ncol = 1))



#################################  End of File  #################################
