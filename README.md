
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ppmm

<!-- badges: start -->
<!-- badges: end -->

The ppmm package implements the proxy pattern-mixture model (PPMM) for
means, proportions, and linear regression coefficients as described in:

(means) Andridge RR, Little RJA (2011). Proxy pattern-mixture analysis
for survey nonresponse. Journal of Official Statistics, 27; 153-180.

(proportions) Andridge RR, Little RJA. (2020) Proxy pattern-mixture
analysis for a binary variable subject to nonresponse. Journal of
Official Statistics, 36; 703-728.

(regression coefficients) West BT, Little RJ, Andridge RR, Boonstra P,
Ware EB, Pandit A, Alvarado-Leiton F (2021). Assessing selection bias in
regression coefficients estimated from nonprobability samples with
applications to genetics and demographic surveys. Annals of Applied
Statistics, 15; 1556-1581.

## Installation

You can install the development version of ppmm from
[GitHub](https://github.com/) with:

``` r
# install package `ppmm`
devtools::install_github("randridge/ppmm")
library(ppmm)
```

## Example

NEEDS MODIFYING

This is a basic example which shows you how to solve a common problem:

``` r
library(ppmm)
#> Warning in .recacheSubclasses(def@className, def, env): undefined subclass
#> "ndiMatrix" of class "replValueSp"; definition not updated
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.
