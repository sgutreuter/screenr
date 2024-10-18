################################################################################
##     R Script:  Sweave_it.R
##
##  Description:  Weave screenr_Tutorial.Rnw
################################################################################
library(tools)
workpath <- file.path(Sys.getenv("DEVEL"), "screenr_maintenance/vignettes")
setwd(workpath)
Sweave("screenr_Tutorial.Rnw")
