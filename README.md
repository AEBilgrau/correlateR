correlateR
==========
[![Build Status](https://api.travis-ci.org/AEBilgrau/correlateR.svg?branch=master)](https://travis-ci.org/AEBilgrau/correlateR)

The R-package `correlateR` features fast, robust, and efficient (as well as inefficient) marginal, partial, semi-partial correlations and covariances of arbitrary conditional order. A good discussion and explanation of marginal (unconditioned), partial, and semi-partial (or, part) correlations is found [here.](http://luna.cas.usf.edu/~mbrannic/files/regression/Partial.html) Another good resource is found [here.](http://www.johndcook.com/blog/2008/11/05/how-to-calculate-pearson-correlation-accurately/)

The package is designed to perform well in both high and low dimensional cases as well as  both on dense and sparse matrices.

The packages is features (or, is planned to feature):
* Covariance/Correlation ✓
* Cross-covariance/Cross-correlation ✓
* Sparse estimation of covariance/correlation
* Sparse estimation of cross-covariance/correlation
* Arbitrary order conditional/partial covariance/correlation
* Arbitrary order conditional/partial cross-covariance/correlation (?)
* Part (or semi-partial) covariances and correlation (?)
* Sparse estimation methods
* Robust estimation methods
* Interface using formulas!
* ... and more!


Naming conventions and interface
--------------------------------
To easily navigate the package some naming conventions has been decided upon.

Lower-case `x`, `y`, `z` always denotes `numeric` vectors while the upper-case counterparts `X`, `Y`, `Z` denote `numeric` matrices where observations correspond to rows and variables/feature to columns. The `Z` and `z` always express the variables conditioned on. Furthermore, `S` is used to denote the empirical (marginal) covariance matrix.

Function names are in camelCase except for some special cases. Otherwise `cor` is for correlation `cov` is for covariance. These are prefixed with `x` or `p` (or both) to denote cross or partial correlations/covariance respectively. For example, `pcor` is the partial correlation and `pxcov` is the partial cross covariance. The available base functions are `cor`, `xcor`, `pcor`, `xpcor`, `cov`, `xcov`, `pcov`, `xpcov`.


Installation
------------
If you wish to install the latest version of `correlateR` directly from the master branch here at GitHub, run 

```R
#install.packages("devtools")  # Uncomment if devtools is not installed
devtools::install_github("AEBilgrau/correlateR")
```

The package is still under heavy development and should be considered unstable. Be sure that you have the [package development prerequisites](http://www.rstudio.com/ide/docs/packages/prerequisites) if you wish to install the package from the source.

*NOTE* The interface and function names may still see significant changes and
modifications!


Alternatives
------------
Some alterative packages are: 
* `corpcor`: Only features estimation of the full partial correlations.
* `ppcor`: 

------------
