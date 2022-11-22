
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FAIMSpredictor

<!-- badges: start -->
<!-- badges: end -->

Let’s use FAIMSpredictor to predict peptide’s most suitable CV settings.

## Installation

You can install the development version of FAIMSpredictor like so:

``` r
devtools::install_github("weixiandeng/FAIMSpredictor")
```
## Download the two machine learning models to your working directory.
## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(FAIMSpredictor)

CV_calc_no_empirical(seq="AAAAAAAAAAAAA",charge=2)
```

In which, seq argument takes peptide sequence as a string, and charge argument takes charge state of the peptide as a numeric variable.
``` r
CV_calc_with_empirical(seq="AAAAAAAAAAAAA",charge=2, cv=45)
```
Similarly, in the second function, once we have a experimentally measured CV value where the maximum intensity of a peptide is observed, we can feed the value to the third argument.

