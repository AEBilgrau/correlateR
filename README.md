correlateR
==========
[![Build Status](https://api.travis-ci.org/AEBilgrau/correlateR.svg?branch=master)](https://travis-ci.org/AEBilgrau/correlateR)

The R-package `correlateR` features fast, robust, and efficient marginal, partial, semi-partial correlations and covariances of arbitrary conditional order. A good discussion and explanation of marginal (unconditioned), partial, and semi-partial (or, part) correlations is found [here.](http://luna.cas.usf.edu/~mbrannic/files/regression/Partial.html) Another good resource is found [here.](http://www.johndcook.com/blog/2008/11/05/how-to-calculate-pearson-correlation-accurately/)

The package is designed to perform well in both high and low dimensional cases as well as  both on dense and sparse matrices.

The packages is planned to feature:
* covariance/correlation
* cross-covariance/correlation
* arbitrary order conditional covariance/correlation
* (?) arbitrary order conditional cross-covariance/correlation
* sparse estimation methods
* robust estimation methods
* ... and more!

Installation
------------
If you wish to install the latest version of `correlateR` directly from the master branch here at GitHub, run 

```R
#install.packages("devtools")  # Uncomment if devtools is not installed
devtools::install_github("AEBilgrau/correlateR")
```

The package is still under heavy development and should be considered unstable. Be sure that you have the [package development prerequisites](http://www.rstudio.com/ide/docs/packages/prerequisites) if you wish to install the package from the source.


Alternatives
------------
Some alterative packages are: `corpcor`, `ppcor`


---
