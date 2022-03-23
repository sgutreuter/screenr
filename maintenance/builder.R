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
library(knitr)
library(tools)

pkgpath <- file.path(Sys.getenv("DEVEL"), "screenr" )
codepath <- file.path(pkgpath, "R")
workpath <- file.path(pkgpath, "maintenance")
manpath <- file.path(pkgpath, "man")
datapath <- file.path(pkgpath, "data")
setwd(pkgpath)

#################################################################################
## Simulate package loading
#################################################################################
##load_all(pkgpath)

#################################################################################
## Compile the vignette
#################################################################################
demo_ <- file.path(pkgpath, "vignettes/screenr_Tutorial.Rmd")
## refs <- read.bib(file.path(pkgpath, "vignettes/screenr_refs.bib")
## keys <- names(refs)
rmarkdown::render(input = demo_)

#################################################################################
## Create the INDEX file
#################################################################################
tools::Rdindex(manpath, outFile = file.path(pkgpath, "INDEX"))



#################################################################################
## Always (re)generate NAMESPACE and documentation, and then check package from
## source tree
#################################################################################
if("package:screenr" %in% search( )) detach("package:screenr")
if(file.exists(file.path(pkgpath, "NAMESPACE"))) {
    file.remove(file.path(pkgpath, "NAMESPACE"))}
use_mit_license()
devtools::document()
devtools::check()
devtools::build_vignettes( )

#################################################################################
## Build README.md
#################################################################################
devtools::build_readme(workpath)
render(file.path(workpath,  "README.Rmd"))
#################################################################################
## Build the package manual (pdf)
#################################################################################
build_manual(pkg = pkgpath, path = pkgpath)

#################################################################################
## Install screenr from the source tree
#################################################################################
if("package:screenr" %in% search( )) detach("package:screenr")
install.packages(pkgpath, repos = NULL, type = "source")

#################################################################################
## Install screenr from GitHub
#################################################################################
if("package:screenr" %in% search( )) detach("package:screenr")
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
