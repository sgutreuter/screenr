################################################################################
##     R Script: install_screenr.R
##
##      Package: screenr
##
##  Description: Install screenr various ways
##
##       Author: Steve Gutreuter
##               E-mail:  sgutreuter@gmail.gov
################################################################################
library(devtools)

pkgpath <- file.path(Sys.getenv("DEVEL"), "screenr" )
codepath <- file.path(pkgpath, "R")
workpath <- file.path(pkgpath, "maintenance")
manpath <- file.path(pkgpath, "man")
datapath <- file.path(pkgpath, "data")
setwd(pkgpath)

################################################################################
## Install from local package source
################################################################################
if("package:screenr" %in% search( )) detach("package:screenr")
install.packages(pkgpath, repos = NULL, type = "source")

################################################################################
## Install from GitHub
################################################################################
## Branch: master
if("package:screenr" %in% search( )) detach("package:screenr")
install_github("sgutreuter/screenr", ref = "master", force = TRUE,
               build_vignettes = TRUE, build_manual = TRUE)

## Branch: v0.5-working
if("package:screenr" %in% search( )) detach("package:screenr")
install_github("sgutreuter/screenr", ref = "v0.5-working", force = TRUE)
