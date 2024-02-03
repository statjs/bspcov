# bspcov

![license](https://img.shields.io/badge/Licence-GPL--2-blue.svg)

An R package for Bayesian Sparse Estimation of a Covariance Matrix

![S&P 500 Example](./figs/thresPPPheatmap.png?raw=true "SP 500 Example")

## Building

To build the package from source, you need to have the following:

```R
# lock the renv
pkgs <- c("...")
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

You can install the development version of `bspcov` from [GitHub](https://github.com/statjs/bspcov) with:

```r
# install.packages("devtools")
devtools::install_github("statjs/bspcov", ref = "main")
```

## Related publications

Lee, Jo, and Lee (2022). The beta-mixture shrinkage prior for sparse covariances with posterior near-minimax rate, Journal of Multivariate Analysis, 192, 105067.  
Lee, Jo, and Lee (2023+). Scalable and optimal Bayesian inference for sparse covariance matrices via screened beta-mixture prior.  
Lee, Lee, and Lee (2023+). Post-processes posteriors for banded covariances, Bayesian Analysis, DOI: 10.1214/22-BA1333.  
Lee and Lee (2023). Post-processed posteriors for sparse covariances, Journal of Econometrics, 236(3), 105475.

## Acknowledgement

This work was supported by the National Research Foundation of Korea(NRF) grant funded by the Korea government(MSIT)   
(RS-2023-00211979, NRF-2022R1A5A7033499, NRF-2020R1A4A1018207, and NRF-2020R1C1C1A01013338)
