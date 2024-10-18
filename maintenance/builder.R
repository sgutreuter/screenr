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
library(qpdf)

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
## Always (re)generate NAMESPACE and documentation, and then check package from
## source tree
#################################################################################
if("package:screenr" %in% search( )) detach("package:screenr")
if(file.exists(file.path(pkgpath, "NAMESPACE"))) {
    file.remove(file.path(pkgpath, "NAMESPACE"))}
use_mit_license()
devtools::document()
tools::Rdindex(manpath, outFile = file.path(pkgpath, "INDEX"))
devtools::check()
##devtools::check(cran = FALSE)
devtools::check(vignettes = FALSE )
devtools::build_vignettes( )
#################################################################################
## NOTE: To build the PDF version of the vignette, see
##       ~/Devel/screenr_maintenance/vignettes/screenr_Tutorial.Rnw.
##       Yes, the duplication is awkward and Rmarkdown is deficient!
#################################################################################
#################################################################################
## Build README.md
#################################################################################
##devtools::build_readme(file.path(pkgpath, "inst/README.Rmd"), path = pkgpath)
##render(file.path(pkgpath, "inst/README.Rmd"), output_file = file.path(pkgpath, "README.md"))
#################################################################################
## Build the package manual (pdf)
#################################################################################
build_manual(pkg = pkgpath, path = pkgpath)

#################################################################################
## Install screenr from the source tree
#################################################################################
if("package:screenr" %in% search( )) detach("package:screenr")
##install.packages(pkgpath, repos = NULL, type = "source")
install_local(pkgpath, repos =  NULL, build_vignettes =  TRUE,
              build_manual =  TRUE, force = TRUE)
#################################################################################
## Install screenr from GitHub
#################################################################################
if("package:screenr" %in% search( )) detach("package:screenr")
devtools::install_github("sgutreuter/screenr", ref = "master", force = TRUE,
                         build_vignettes = TRUE, build_manual = TRUE)
#################################################################################
## Reload screenr (if needed)
#################################################################################
reload("screenr")

#################################################################################
## Remove screenr package
#################################################################################
##remove.packages("screenr")
