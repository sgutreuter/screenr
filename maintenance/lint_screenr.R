#################################################################################
##      R PROGRAM: lint_screenr.R
##
##        PROJECT: screenr Package
##
##    DESCRIPTION: Lint the entire screenr package
##
##     WRITTEN BY: Steve Gutreuter
##                 E-mail:  sgutreuter@gmail.gov
#################################################################################

#################################################################################
## Set paths and working directory
codepath <- file.path(Sys.getenv("DEVEL"), "screenr/R")
workpath <- file.path(Sys.getenv("DEVEL"), "screenr/maintenance")
setwd(workpath)

#################################################################################
## Lint the screenr package
#################################################################################
library(lintr)
lint_dir(codepath)

#################################  End of File  #################################
