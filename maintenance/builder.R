################################################################################
##     R Script: builder.R
##
##      Package: screenr
##
##  Description: Code for building the screenr package
##
##       Author: Steve Gutreuter
##               E-mail:  sgutreuter@gmail.gov
################################################################################
rm(list =  ls())

library(devtools)
library(roxygen2)
library(rmarkdown)

pkgpath <- file.path(Sys.getenv("DEVEL"), "screenr" )
codepath <- file.path(pkgpath, "R")
workpath <- file.path(pkgpath, "maintenance")
datapath <- file.path(pkgpath, "data")
setwd(pkgpath)

#################################################################################
## Simulate package loading
#################################################################################
##load_all(pkgpath)

#################################################################################
## Compile the vignette
#################################################################################
## library(bibtex)
## demo <- file.path(pkgpath, "vignettes/screenr_demo.Rmd")
## refs <- read.bib(file.path(pkgpath, "vignettes/screenr_refs.bib")
## keys <- names(refs)
## rmarkdown::render(input = demo)


#################################################################################
## Always (re)generate NAMESPACE and documentation, and then check package from
## source tree
#################################################################################
if("package:screenr" %in% search( )) detach("package:screenr")
if(file.exists(file.path(pkgpath, "NAMESPACE"))) {
    file.remove(file.path(pkgpath, "NAMESPACE"))}
devtools::document()
devtools::check()

#################################################################################
## Build vignette
#################################################################################
devtools::build_vignettes( )

#################################################################################
## Build the package manual (pdf)
#################################################################################
##build_manual(as.package(pkgpath))   ## Fails:  debugonce(build_manual)


#################################################################################
## Install screenr from the source tree
#################################################################################
install.packages(pkgpath, repos = NULL, type = "source")

#################################################################################
## Install screenr from GitHub
#################################################################################
devtools::install_github("sgutreuter/screenr", ref = "master")
devtools::install_github("sgutreuter/screenr", ref = "testing")

#################################################################################
## Reload screenr (if needed)
#################################################################################
reload("screenr")

#################################################################################
## Remove screenr package
#################################################################################
##remove.packages("screenr")
