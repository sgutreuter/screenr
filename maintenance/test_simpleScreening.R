#################################################################################
##      R PROGRAM: simpleScreening_testing.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: Testing sandbox for simpleScreening
##
##     WRITTEN BY: Steve Gutreuter
##                 E-mail:  sgutreuter@gmail.gov
#################################################################################

#################################################################################
## Set paths and working directory
#################################################################################
codepath <- file.path(Sys.getenv("DEVEL"), "screenr/R")
workpath <- file.path(Sys.getenv("DEVEL"), "screenr/maintenance")
setwd(workpath)

library(pROC)
library(tidyverse)

#################################################################################
## Source the screenr R code for experimentation and testing
#################################################################################
source(file.path(codepath, "helperFunctions.R"))
source(file.path(codepath, "simpleScreening.R"))
source(file.path(workpath, "genData.R"))
#################################################################################
## or attach the package
#################################################################################
library(screenr)

#################################################################################
## Test code
#################################################################################
## simpleScreening
simpleScreening <- function(formula, data){
    mf <- match.call(expand.dots = FALSE)
    call <- mf
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    dat <- data[complete.cases(data[, names(mf), ]), ]
    y <- model.response(mf, "numeric")
    if(!all(y %in% c(0, 1))) stop("Response variable must be binary (0, 1)")
    prev <- mean(y, na.rm = TRUE)
    preds <- mf[, -1]
    npreds <- dim(preds)[2]
    score <- apply(preds, 1, sum)
    is.roc <- pROC::roc(y, score, auc = TRUE, direction = "<")
    scores <- cbind(dat, score = score)
    result <- list(Call = call,
                   Prevalence = prev,
                   ISroc = is.roc,
                   Scores = scores)
    class(result) <- "simplescreenr"
    invisible(result)
}

## summary.simplescreenr
summary.simplescreenr <- function(object, ...){
    if(!("simplescreenr" %in% class(object)))
        stop("object not a simplescreenr object")
    cat("Call:\n")
    print(object$Call)
    cat("\nPrevalence (In-sample prevalence of condition):\n")
    print(object$Prevalence)
    cat("\nReceiver Operating Characteristics:\n")
    auc <- round(as.numeric(object$ISroc$auc), digits = 4)
    cat(paste("\nIn-sample area under the ROC curve: ",
              auc, "\n", sep = ""))
}



## plot.simplescreenr
plot.simplescreenr <- function(x, plot_ci = TRUE, print_ci = TRUE,
                               conf_level = 0.95, bootreps = 2000,...){
    if(!class(x) == "simplescreenr") stop("x is not a simplescreenr object")
    plt <- plot(x$ISroc, print.auc = TRUE, ...)
    if(plot_ci | print_ci){
        ciplt <- pROC::ci.thresholds(x$ISroc, boot.n = bootreps,
                                     progress = "none",
                                     conf.level = conf_level,
                                     thresholds = "local maximas")
        }
    if(print_ci){
        threshold <- as.numeric(rownames(ciplt$sensitivity)) + 0.5
        citable <- data.frame(cbind(threshold, ciplt$sensitivity,
                                    ciplt$specificity))
        names(citable) <- c("threshold", "se.low", "se.median",
                            "se.high", "sp.low", "sp.median", "sp.high")
        row.names(citable)  <- 1:length(threshold)
    }
    if(plot_ci) plot(ciplt)
    if(print_ci) return(citable)
}


## print.simplescreenr
print.simplescreenr <- function(x, quote = FALSE, ...){
    if(!("simplescreenr" %in% class(x))) stop("x not a simplescreenr object")
    cat("\nIn-sample (overly optimistic) sensitivity and specificity:\n")
    df_ <- pROC::coords(x$ISroc, transpose = FALSE)
    df_["threshold"] <- df_["threshold"] + 0.5
    print(df_)
}

## getROC
getROC <- function(x, simplify = TRUE){
    sensitivity <- specificity <- NULL
    if(!class(x) %in% c("binomscreenr", "simplescreenr"))
        stop("x not a binomscreenr or simplescreenr object.")
    if(class(x) == "binomscreenr"){
        res <- pROC::coords(x$CVroc, transpose = FALSE)
    } else {
        res <- pROC::coords(x$ISroc, transpose = FALSE)
        res$threshold  <-  res$threshold + 0.5
    }
    if(simplify) {
        cleaned <- res %>%
            dplyr::group_by(sensitivity) %>%
            dplyr::summarize(specificity = max(specificity))
        res <- dplyr::right_join(res, cleaned)
    }
    res
}

################################################################################
## Testing
################################################################################
smpl <- simpleScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5 ,
                        data = unicorns)

str(smpl)

summary(smpl)

print(smpl)

(wtf <- getROC(smpl))
(wtf <- getROC(smpl, simplify = FALSE))

debugonce(plot.simplescreenr)
(ci <- plot(smpl))
plot(smpl, print_ci = FALSE, type = "S")



#################################  End of File  #################################
