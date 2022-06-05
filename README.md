
<!-- README.md is generated from README.Rmd -->

# screenr

**An R package to enable screening in/out subjects who are likely to to
test positive/negative, respectively.**

## Description

Package screenr enables easy development and validation of diagnostic
test screening tools. It is designed to enable those with only a basic
familiarity with R to develop, validate and screening tools for
diagnostic tests.

Consider the situation where a definitive test for some condition is
expensive, and the condition is rare. In that case, universal testing
would not be efficient in terms of the yield of postive results per test
performed. Now suppose that responses to a set of simple diagnostic
questions or observations may be predictive of the definitive test
result. Package screenr enables the estimation of thresholds for making
decisions about when to perform the definitive test on newly observed
subjects based on Receiver Operating Characteristics (ROC) estimated
from an initial sample. The choice of a particular screening threshold
is left to the user, and should be based on careful consideration of
application-specific tradeoffs between sensitivity (true positive
fraction) and specificity (true negative fraction).

## Installation

### Installing from GitHub in an R session

Open an R session and do:

``` r
## Install devtools if not already installed:
install.packages("devtools")
## Install screenr:
devtools::install_github("sgutreuter/screenr", build_vignettes = TRUE,
                         build_manual = TRUE)
```

**or**

### Download and install from a command line

1.  Browse to <https://github.com/sgutreuter/screenr>.

2.  Click on the **Code** pull-button near the upper-right corner,
    select **Download ZIP** and save the zip archive in a location of
    your choice.

3.  Unzip the archive into an empty temporary directory.

4.  Open a shell/terminal in the directory containing the unzipped
    screenr directory.

5.  Execute ‘R CMD INSTALL ./screenr’

## Getting started

The tutorial is recommended for first-time users:

``` r
    library(screenr)
    vignette("screenr_Tutorial", package = "screenr")
```

## Bug reports

<https://github.com/sgutreuter/screenr/issues>
