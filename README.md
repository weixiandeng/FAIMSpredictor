
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FAIMSpredictor

<!-- badges: start -->
<!-- badges: end -->

Let’s use FAIMSpredictor to predict peptide’s most suitable CV settings.

## Installation

You can install the development version of FAIMSpredictor like so:

``` r
devtools::install_github("weixiandeng/FAIMS_CV_predictor")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(FAIMSpredictor)

CV_calc_no_empirical("AAAAAAAAAAAAA",2)

CV_calc_with_empirical("AAAAAAAAAAAAA",2, 45)
```
