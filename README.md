# bspcov

[![license](https://img.shields.io/badge/License-GPL--2-blue.svg)](https://github.com/statjs/bspcov/blob/main/LICENSE.txt)
[![R](https://img.shields.io/badge/R-%3E%3D4.2-blue.svg)](https://cran.r-project.org/package=bspcov)
[![release](https://img.shields.io/github/v/release/statjs/bspcov)](https://github.com/statjs/bspcov/releases/latest)
[![Check R packages](https://github.com/statjs/bspcov/actions/workflows/check.yml/badge.svg)](https://github.com/statjs/bspcov/actions/workflows/check.yml)
[![Deploy R packages](https://github.com/statjs/bspcov/actions/workflows/deploy.yml/badge.svg)](https://github.com/statjs/bspcov/actions/workflows/deploy.yml)  
[![dev-latest](https://img.shields.io/github/release-date-pre/statjs/bspcov?filter=dev-latest&label=dev-latest)](https://github.com/statjs/bspcov/releases/tag/dev-latest)
[![Check R packages (dev)](https://github.com/statjs/bspcov/actions/workflows/dev-check.yml/badge.svg)](https://github.com/statjs/bspcov/actions/workflows/dev-check.yml)
[![Deploy R packages (dev)](https://github.com/statjs/bspcov/actions/workflows/dev-deploy.yml/badge.svg)](https://github.com/statjs/bspcov/actions/workflows/dev-deploy.yml)

An R package for Bayesian Sparse Estimation of a Covariance Matrix

![S&P 500 Example](https://github.com/statjs/bspcov/blob/main/figs/thresPPPheatmap.png?raw=true "SP 500 Example")

## Building

To build the package from source, you need to have the following:

```R
# lock the renv
pkgs <- c("GIGrvg", "coda", "progress", "BayesFactor", "MASS", "mvnfast", 
          "matrixcalc", "matrixStats", "purrr", "dplyr", "RSpectra", "Matrix", 
          "plyr", "CholWishart", "magrittr", "future", "furrr", "ks", "ggplot2", 
          "ggmcmc", "caret", "FinCovRegularization", "mvtnorm", "stats", 
          "patchwork", "reshape2", "future.apply")
renv::snapshot(packages = pkgs)

# update docs
devtools::document()
```

```bash
## check package
VERSION=$(git describe --tags | sed 's/v//g')

## build manual
R CMD Rd2pdf --force --no-preview -o bspcov-manual.pdf .

## build package
sed -i '' "s/Version: [^\"]*/Version: ${VERSION}/g" "DESCRIPTION"
R CMD build .
```

## Installation

You can install the `bspcov` package from [CRAN](https://CRAN.R-project.org/package=bspcov):
```r
install.packages("bspcov")
```
or the development version from GitHub, by using the function `install_github()` from [`devtools`](https://CRAN.R-project.org/package=devtools) package:
```r
devtools::install_github("statjs/bspcov", ref = "main")
```

## Related publications

* Lee, K., Jo, S., & Lee, J. (2022). The beta-mixture shrinkage prior for sparse covariances with near-minimax posterior convergence rate. Journal of Multivariate Analysis, 192, 105067, DOI: 10.1016/j.jmva.2022.105067.
* Lee, K., Jo, S., Lee, K., & Lee, J. (2024). Scalable and optimal Bayesian inference for sparse covariance matrices via screened beta-mixture prior. Bayesian Analysis, 1(1), 1-28, DOI: 10.1214/24-BA1495.
* Lee, K., Lee, K., & Lee, J. (2023). Post-processed posteriors for banded covariances. Bayesian Analysis, 18(3), 1017-1040, DOI: 10.1214/22-BA1333.
* Lee, K., & Lee, J. (2023). Post-processed posteriors for sparse covariances. Journal of Econometrics, 236(1), 105475, DOI: 10.1016/j.jeconom.2023.105475.

## Acknowledgement

This work was supported by the National Research Foundation of Korea(NRF) grant funded by the Korea government(MSIT)   
(RS-2023-00211979, NRF-2022R1A5A7033499, NRF-2020R1A4A1018207, and NRF-2020R1C1C1A01013338)
