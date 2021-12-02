################################################################################
##      R PROGRAM: builder.R
##
##    DESCRIPTION: Code for building the screenr package
##
##     WRITTEN BY: Steve Gutreuter
##                 E-mail:  sgutreuter@gmail.gov
################################################################################
library(devtools)
library(roxygen2)

pkgpath <- file.path(Sys.getenv("DEVEL"), "screenr" )
codepath <- file.path(pkgpath, "R")
workpath <- file.path(pkgpath, "maintenance")
datapath <- file.path(pkgpath, "data")
setwd(pkgpath)

#################################################################################
## Always (re)generate documentation, and then check package from source tree
#################################################################################
devtools::document()
devtools::check()

#################################################################################
## Install screenr from the source tree
#################################################################################
install.packages(workpath, repos = NULL, type = "source")

#################################################################################
## Install screenr from GitHub
#################################################################################
devtools::install_github("sgutreuter/screenr", ref = "master")
devtools::install_github("sgutreuter/screenr", ref = "testing")

#################################################################################
## Remove screenr package
#################################################################################
remove.packages("screenr")
