#################################################################################
##      R PROGRAM: test_screenr_pkg.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: Testing sandbox code for screenr
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
## Either source the screenr R code or load the package, but not both:
#################################################################################
## 1. Source the code
## source(file.path(codepath, "binomialScreening.R"))
## source(file.path(codepath, "geebinomScreening.R"))
## source(file.path(codepath, "helperFunctions.R"))
## source(file.path(codepath, "mebinomScreening.R"))
## source(file.path(codepath, "simpleScreening.R"))
## source(file.path(workpath, "genData.R"))

##   -- or --

## 2. Load the screener package
 library(screenr)
## ?screenr
 data(screenr::unicorns)
## help(unicorns)

#################################################################################
## Create new data for prediction
#################################################################################
new <- data.frame(ID = c('"Bernie P."', '"Alice D."'), Q1 = c(0, 0), Q2 = c(0, 0),
                  Q3 = c(1, 0), Q4 = c(0, 0), Q5 = c(1, 0))


#################################################################################
## binomScreening
#################################################################################

?binomialScreening
bsobj <- binomialScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5,
                           data = unicorns, link = "logit", Nfolds = 20)
summary(bsobj)
plot(bsobj)
print(bsobj)
?testCounts
testCounts(bsobj)
getROC(bsobj)
## Compute point estimates of their predicted probabilities testing positive:
predict(bsobj$ModelFit, newdata = new, type = "response")
## or, compute directly
?inverseLink
inverseLink("logit",
            as.matrix(cbind(rep(1, nrow(new)), new[, 2:6])) %*%
                            as.matrix(bsobj$ParamEst, ncol = 1))

#################################################################################
## geebinomScreening
#################################################################################
geeobj <- geebinomScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5, id = "clinic",
                            data = unicorns, link = "logit")
## Error in eval(extras, data, env) : object 'id' not found
summary(geeobj)
plot(geeobj)
## But the following works:
wtf <-  gee::gee(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5, id = clinic,
                        data = unicorns, family = binomial(link = "logit"))


geewrapper1 <- function(formula, id, data, family){
    ida <- substitute(id)
    f <- parent.frame()
    result <- geepack::geeglm(formula, id = eval(ida, envir = data),
                              data = data, family = family)
    result
}
debugonce(geewrapper1)
geewrapper1(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5, id = clinic, data = unicorns,
                        family = binomial(link = "logit"))

geewrapper2 <- function(formula, id, data, family){
    family <- substitute(family)
    id <- substitute(id)
    result <- gee::gee(formula, id = id, data = data, family = family)
    result
}
debugonce(geewrapper2)
geewrapper2(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5, id = clinic, data = unicorns,
            family = binomial(link = "logit"))

str(geeobj)

summary(geeobj)

print(geeobj)

getROC(geeobj)
getROC(gee, simplify = FALSE)

plot(geeobj)
plot(geeobj, print_ci = FALSE)


#################################################################################
## mebinomScreening
#################################################################################
?mebinomScreening
meobj <- mebinomScreening(testresult ~ Q1 + Q2 + Q3 + Q4 + Q5, id = "clinic",
                          data = unicorns, link = "logit")

str(meobj)

summary(meobj)

print(meobj)

getROC(meobj)
getROC(meobj, simplify = FALSE)

plot(meobj)
plot(meobj, print_ci = FALSE)

#################################################################################
## simpleScreening (The last resort.  Ugghhh!!!)
#################################################################################

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
