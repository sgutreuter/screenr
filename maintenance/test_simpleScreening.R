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
codepath <- file.path(Sys.getenv("DEVEL"), "screenr/R")
workpath <- file.path(Sys.getenv("DEVEL"), "screenr-testing")
setwd(workpath)

#################################################################################
## Source the screenr R code for experimentation and testing
#################################################################################
source(file.path(codepath, "helperFunctions.R"))
source(file.path(codepath, "simpleScreening.R"))
source(file.path(workpath, "genData.R"))

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
    is.roc <- pROC::roc(y, score, auc = TRUE)
    scores <- cbind(dat, score = score)
    result <- list(Call = call,
                   Prevalence = prev,
                   ISroc = is.roc,
                   Scores = scores)
    class(result) <- "simplescreenr"
    invisible(result)
}

## plot.simplescreenr
plot.simplescreenr <- function(x, plot_ci = TRUE, print_ci = TRUE,
                               conf_level = 0.95, bootreps = 2000,...){
    if(!class(x) == "simplescreenr") stop("x is not a simplescreenr object")
    plt <- plot(x$ISroc, print.auc = TRUE, ...)
    if(plot_ci | print_ci){
        ciplt <- ci.thresholds(x$ISroc, boot.n = bootreps, progress = "none",
                               conf.level = conf_level,
                               thresholds = "local maximas")
        }
    if(print_ci){
        threshold <- 0:(length(dimnames(ciplt$sensitivity)[[1]]) -1 )
        citable <- data.frame(cbind(threshold, ciplt$sensitivity,
                                    ciplt$specificity))
        names(citable) <- c("threshold", "se.low", "se.median",
                            "se.high", "sp.low", "sp.median", "sp.high")
        row.names(citable)  <- threshold + 1
    }
    if(plot_ci) plot(ciplt)
    if(print_ci) return(citable)
}

## print.simplescreenr
print.simplescreenr <- function(x, quote = FALSE, ...){
    if(!("simplescreenr" %in% class(x))) stop("x not a simplescreenr object")
    cat("\nIn-sample (overly optimistic) sensitivity and specificity:\n")
    df_ <- data.frame(score = 0:(length(x$ISroc$sensitivities) - 1),
                      sensitivity = x$ISroc$sensitivities,
                      specificity = x$ISroc$specificities)
    print(df_)
}
print(smpl, quote = TRUE)

## getROC
getROC <- function(x, simplify = TRUE){
    sensitivity <- specificity <- NULL
    if(!class(x) %in% c("binomscreenr", "simplescreenr"))
        stop("x not a binomscreenr or simplescreenr object.")
    if(class(x) == "binomscreenr"){
        obj <- x$CVroc
        th <- obj$thresholds
    } else {
        obj <- x$ISroc
        th <- 0:(length(x$ISroc$sensitivities) - 1)
    }
    res <- data.frame(threshold = th,
                      sensitivity = obj$sensitivities,
                      specificity = obj$specificities)
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

getROC(smpl)
getROC(smpl, simplify = FALSE)

plot(smpl)
plot(smpl, print_ci = FALSE)



#################################  End of File  #################################
