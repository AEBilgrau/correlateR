correlateR
==========
[![Build Status](https://api.travis-ci.org/AEBilgrau/correlateR.svg?branch=master)](https://travis-ci.org/AEBilgrau/correlateR)

The R-package `correlateR` features fast, robust, and efficient (as well as inefficient) marginal, partial, semi-partial correlations and covariances of arbitrary conditional order. A good discussion and explanation of marginal (unconditioned), partial, and semi-partial (or, part) correlations can be found [here.](http://luna.cas.usf.edu/~mbrannic/files/regression/Partial.html) Another good resource is found [here.](http://www.johndcook.com/blog/2008/11/05/how-to-calculate-pearson-correlation-accurately/)

The package is designed to perform well in both high and low dimensional cases as well as both on dense and sparse matrices.

The packages is features (or, is planned to feature):

 - [x] `cor`/`cov` Marginal (unconditional) correlation/covariance. These basic 
       functions can be prefixed to yield other correlation/covariance 
       estimates. This covariance is also known as the auto-correlation, the 
       variance-covariance, or simply the variance (in the generalized sense).
    - [x] `p`-prefix: *p*artial (arbitrary order) correlation and covariance.
    - [x] `x`-prefix: *cross* correlation and covariance.
    - [ ] `P`-prefix: *P*art (semi-partial) correlation and covariances (?)
    - [ ] `s`-prefix: *s*parse shrinkage estimation methods (?)
    - [ ] `r`-prefix: *r*obust estimation methods. E.g. Minimum Covariance 
          Determinant
    - [x] `S`-prefix: *S*hrinkage estimation. (Or, `d` for *d*ense shrinkage?) 
 - [ ] Interface using formulas `~`.
 - [ ] Conversion between `cov` and `cor` and `pcor` functions.
    - [ ] `cov2cor`, `cor2cov`, `cor2pcor`, `pcor2cor`(?)
 - [ ] Conditional and unconditional independence test
    - [ ] `cor.text`, `pcor.test`
    - [ ] Also with cross, sparse, shrinked, robust, etc., versions
 - [ ] Canonical correlation analysis (CCA)
    - [ ] Also with cross, sparse, shrinked, robust, etc., versions
 - [ ] `pre` (alternative to `cov`) direct estimation of the precision matrix
       or concentration matrix.
    - [ ]  Also with cross, sparse, robust, etc., versions
 - [ ] ... and more! (??)
 
Hence the following core-functons are available:
 - [x] `xcor` Cross-correlation
 - [x] `xcov` Cross-covariance
 - [x] `pcor` Partial correlation (arbitrary order)
 - [x] `pcov` Partial covariance (arbitrary order)
 - [x] `pxcor` Partial cross-correlation (arbitrary order)
 - [x] `pxcov` Partial cross-covariance (arbitrary order)
 - [ ] `scor` Sparse correlation
 - [ ] `scov` Sparse covariance
 - [ ] `sxcor` Sparse cross-correlation
 - [ ] `sxcov` Sparse cross-covariance
 - [ ] `spcor` Sparse partial correlation (arbitrary order)
 - [ ] `spcov` Sparse partial covariance (arbitrary order)
 - [ ] `spxcor` Sparse partial cross-correlation (arbitrary order)
 - [ ] `spxcov` Sparse partial cross-covariance (arbitrary order)



Naming conventions and interface
--------------------------------
To easily navigate the package some naming conventions has been decided upon.

Lower-case `x`, `y`, `z` always denotes `numeric` vectors while the upper-case counterparts `X`, `Y`, or `Z` denote a `numeric` `matrix` where observations correspond to rows and variables/feature to columns. The `Z` and `z` always express the variables conditioned on. Furthermore, `S` is used to denote the empirical (marginal) covariance matrix.

Function names are in camelCase except for some special cases. Otherwise `cor` is for correlation `cov` is for covariance. These are prefixed with `x` or `p` (or both) to denote cross or partial correlations/covariance respectively. For example, `pcor` is the partial correlation and `pxcov` is the partial cross covariance. 


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


Using correlateR
----------------
A small tutorial goes here!


Alternative packages
--------------------
There are some alternative packages on CRAN form which some inspiration have been drawn. 
* `corpcor`: Only features estimation of the full partial correlations.
* `ppcor`: Partial and semi-partial correlations

### Benchmarking
A benchmarking of the functions goes here.


Reference
---------
References goes here!


--------------------------------------------------------------------------------
