#################################################################################
## screenr Package-building sandbox
#################################################################################

library(devtools)
library(roxygen2)
workpath <- file.path(Sys.getenv("DEVEL"), "screenr")
setwd(workpath)

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
