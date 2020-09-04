#################################################################################
##      R PROGRAM: test_geebinomScreening.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: Testing sandbox for geebinomScreening
##
##     WRITTEN BY: Steve Gutreuter
##                 E-mail:  sgutreuter@gmail.gov
#################################################################################

#################################################################################
## Set paths and working directory
#################################################################################
basepath <- file.path(Sys.getenv("DOCS"), "Computing/Devel/screenr")
codepath <- file.path(basepath, "R")
workpath <- file.path(basepath, "maintenance")
setwd(workpath)

#################################################################################
## Source the screenr R code for experimentation and testing
#################################################################################
source(file.path(workpath, "genData.R"))
source(file.path(codepath, "helperFunctions.R"))
##source(file.path(codepath, "geebinomScreening.R"))

#################################################################################
## Create new data for prediction
#################################################################################
new <- data.frame(ID = c('"Bernie P."', '"Alice D."'), Q1 = c(0, 0), Q2 = c(0, 0),
                  Q3 = c(1, 0), Q4 = c(0, 0), Q5 = c(1, 0))

#################################################################################
## geebinomScreening
#################################################################################

geebinomScreening <- function(formula,
                             id = NULL,
                             corstr = "exchangeable",
                             data = NULL,
                             Nfolds = 40L,
                             ...){
    if(!inherits(formula, "formula")) stop("Specify an model formula")
    if(!is.data.frame(data)) stop("Provide a data frame")
    if(!corstr %in% c("exchangeable", "independence", "ar1", "unstructured",
                      "userdefined"))
        stop("Invalid corstr")
    call <- match.call()
    formx <- update(formula, paste("~ . +", id))
    mf <- stats::model.frame(formx, data = data)
    dat <- eval(mf, parent.frame())
    dat <- dat[complete.cases(dat), ]
    dat <- dat[order(dat[[substitute(id)]]), ]
    if(Nfolds > 0.20*dim(dat)[1])
        stop("Nfolds must be < 20% of number of complete observations")
    y <- stats::model.response(dat, "numeric")
    if(!all(y %in% c(0, 1))) stop("Response variable must be binary (0, 1)")
    prev <- mean(y, na.rm = TRUE)
    geeargs <- list(formula = formx,
                    id = as.name(id),
                    corstr = corstr,
                    data = dat,
                    family = binomial(link = "logit"),
                    ...)
    lrfit <- do.call(geepack::geeglm, geeargs)
    is.roc <- pROC::roc(lrfit$y, as.vector(lrfit$fitted.values))
    N <- nrow(dat)
    holdouts <- split(sample(1:N), 1:Nfolds)
    cv.results <- data.frame(NULL)
    for(i in 1:Nfolds){
        geeargs[["data"]] <- dat[-holdouts[[i]], ]
        res <- do.call(geepack::geeglm, geeargs)
        pred.prob <- inverseLink(predict(res, newdata = dat[holdouts[[i]], ]),
                                 link = "logit")
        y <- model.response(dat[holdouts[[i]], ])
        cv.results <- rbind(cv.results,
                            data.frame(cbind(fold = rep(i, length(pred.prob)),
                                             y = y,
                                             data = dat[holdouts[[i]],],
                                             cv.pred.prob = pred.prob)))
    }
    cv.roc <- pROC::roc(cv.results$y, cv.results$cv.pred.prob,
                                  auc = TRUE, ci = TRUE, of = "sp",
                                  se = seq(0, 1, 0.05), ci.type = "shape")
    class(cv.results) <-  c("cv.predictions", "data.frame")
    result <- list(Call = call,
                   ModelFit = lrfit,
                   Prevalence = prev,
                   ParamEst = lrfit$coefficients,
                   ISroc = is.roc,
                   CVpreds = cv.results,
                   CVroc = cv.roc)
    class(result) <- "binomscreenr"
    invisible(result)
}

################################################################################
## Try it
################################################################################

debugonce(geebinomScreening)
res1 <- geebinomScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5,
                          id = "clinic", data = unicorns, Nfolds = 2L)

summary(res1)

new <- data.frame(ID = c('"Bernie P."', '"Alice D."'), Q1 = c(0, 0),
                   Q2 = c(0, 0), Q3 = c(1, 0), Q4 = c(0, 0), Q5 = c(1, 0),
                   clinic = factor(c("C-5", "C-15"),
                                  levels = c("C-1", "C-10", "C-11", "C-12", "C-13",
                                  "C-14", "C-15", "C-16", "C-17,", "C-18",
                                  "C-19", "C-2", "C-20", "C-3", "C-4",
                                  "C-5", "C-6", "C-7", "C-8", "C-9")))

inverseLink(as.matrix(cbind(rep(1, nrow(new)), new[, 2:6])) %*%
            as.matrix(res1$ParamEst, ncol = 1),
            link = "logit")

predict(res1$ModelFit, newdata = new, type = "response", re.form = ~0)
