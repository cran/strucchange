## strucchange 1.6-0

* Code repository changed from R-Forge to Codeberg at:
  <https://codeberg.org/zeileis/strucchange/>

* Added altdoc page with overview and documentation at:
  <https://zeileis.codeberg.page/strucchange/>

* A new hands-on introduction to the package was written for the
  altdoc page which is also available as:
  `vignette("strucchange-seatbelt", package = "strucchange")`

* Fixed long-standing issue of multiple definitions of `computeEmpProc()`
  within `mefp()` with different arguments. Now all definitions have all
  of the arguments (set to `NULL` by default) and check that only the
  required arguments are specified.

* Similarly the issue of multiple definitions of `plotProcess()`
  with different arguments in `efpFunctional()` has been fixed. Now all
  arguments are the same and differences in the defaults are resolved
  within the different functions, also catching potentially erroneous
  specifications.
  
* Updated `structure()` calls to use `names = ...` instead of `.Names = ...` etc.


## strucchange 1.5-4

* Updated reference output from examples for R CMD check.


## strucchange 1.5-3

* Fixups for R CMD check on different platforms (thanks to support from
  Kurt Hornik and Tomas Kalibera).

* Fixed auto-detection of problems with sorting in `summary()` methods of
  `breakpointsfull` objects (reported by Spencer Graves). In some extreme
  cases with many breakpoints and very short segments, the `sort = TRUE`
  approach of displaying the `breakpoints` fails which was not reported
  correctly automatically.


## strucchange 1.5-2

* Added `lrvar = FALSE` argument to `efp()` so that optionally a long-run variance
  estimator (Andrews or Newey-West) instead of the standard OLS estimator can
  be used for the error variance.

* Added `vcov. = NULL` argument to `efp()` so that optionally other covariance
  matrix estimators can be plugged into RE and ME tests (e.g., `sandwich()`,
  `kernHAC()`, `NeweyWest()`, etc.).

* The default `recresid()` gained an argument `qr.tol = 1e-7` which allows to
  pass another tolerance to `lm.fit()` for detecting linear dependencies in
  (small) subsamples. Furthermore, an argument `engine = c("R", "C")` has been
  added along with an alternative faster C implementation (by Nikolaus
  Umlauf).
  
* The formula method of `breakpoints()` now passes `...` to `recresid()`, e.g., for
  the `qr.tol` and engine arguments above.

* The `breaks` argument of the `breakpoints()` formula argument is now checked
  to be at least 1.

* `breakpoints(..., hpc = "foreach")` now also works if the `foreach` package
  is not attached.

* Improved support for formulas like `y ~ .` in `efp()`, `Fstats()`, and
  `breakpoints()` (suggested by Matthieu Stigler).

* Bug fix in `gefp(..., decorrelate = FALSE)`. Scaling is done with the square
  root of the diagnoal of the variance - as opposed to the diagonal of the
  square root of the variance (reported by Dries Debeer).


## strucchange 1.5-1

* `ordL2BB()` now uses a direct simulation method based on `mvtnorm::rmvnorm()`
  which is much faster, making the computation of p values and critical values
  for the ordinal maxLM statistic much faster and feasible "on the fly".

* Reduced number of significant digits in the `summary()` for `breakpoints` to
  `getOption("digits") - 3`.

* Reference output updated for recent versions of R.


## strucchange 1.5-0

* Added new `efpFunctional` generators for conducting various types of
  structural change tests based on empirical fluctuation processes of class
  `"gefp"`. In particular a (maximum) MOSUM functional was added as well
  as several functionals suitable for aggregation along categorical
  variables. The documentation for previously available functionals
  such as `supLM()` was also enhanced.
  
* The new functionals mentioned above for assessing parameter instability
  along (ordered) categorical variables are `catL2BB` (unordered),
  `ordL2BB` and `ordwmax` (ordered). These are discussed in more detail in
  Merkle, Fan, and Zeileis (2013, Psychometrika).
  
* Added a new default method for `sctest()`. This essentially just calls
  `gefp(object, fit = NULL)` and then (optionally) calls `plot()` and `sctest()`
  using the specified functional. However, several convenience options have
  been added, e.g., using the maximum likelihood information (rather than
  the outer product of gradients) for the covariance matrix or specifying
  the functional via a character string.

* Documentation of the `sctest()` generic and its methods have been enhanced.
  Methods for `formula`, `efp`, and `Fstats` are suitable for assessing
  structural changes in linear regression models while the `default` and
  `gefp` methods (see above) are suitable for general parametric models.

* Improved `plot()` method for `gefp`/`efpFunctional` to allow for more
  flexibility in boundary display. Rather than only `boundary = TRUE` or
  `FALSE` one can now specify a list of graphical parameters, e.g.,
  `boundary = list(col = "slategray", lty = 2)`.

* Updated Depends/Imports in DESCRIPTION/NAMESPACE with new R CMD check
  requirements.


## strucchange 1.4-7

* `plotProcess()` function in `"efpFunctional"` objects now takes a
  `boundary = TRUE` argument by default which can be set to `FALSE` to
  suppress plotting of the boundary function.
  
* Added a check (and a more intelligible warning) in the `formula`
  method of `breakpoints()` whether the `breaks` argument supplied by
  the user is too large.


## strucchange 1.4-6

* Default `recresid()` can now also deal with regressors that do not
  vary across (small) subsamples.


## strucchange 1.4-5

* Further improvements in new `recresid()` default method.
  Now also works correctly if some coefficients are not identified
  on the initial subsamples in the recursion.

* Resaved datasets to reduce storage requirements.

* Fixed bug in `breakpoints()` for time series that contain `NA`s.


## strucchange 1.4-4

* Default `recresid()` method now tries to choose adaptively between
  using the faster updating formula and the slower full
  QR decomposition to yield numerically more stable results.
  In previous versions of the function the QR decomposition was
  used only in the first iteration.

* Improvement in `breakdates()` computations.


## strucchange 1.4-3

* Speed-up in `breakpoints()` for the intercept-only case,
  i.e., `breakpoints(y ~ 1)`.


## strucchange 1.4-2

* Improved time index computations in `gefp()`.

* Added replication notes in `?durab`.


## strucchange 1.4-1

* `efp()`, `Fstats()`, and `breakpoints()` are now more cautios about using
  time series properties from the data and try to check whether any
  `NA`s were removed. In general, the functions will yield best results
  if all `NA` processing is done before calling them.
  
* Better handling of time series properties for the boundaries in
  the examples of `SP2001`.


## strucchange 1.4-0

* Added optional high performance computing support by means of the
  `foreach` package for the `breakpoints()` formula method. This can
  be leveraged to alleviate the computational burden in the dynamic
  programming approach. Simply register a parallel backend (e.g.,
  by means of `doMC` or `doSNOW`) and call `breakpoints()` with
  additional argument `hpc = "foreach"`.


## strucchange 1.3-7

* Added optional `start` end `end` arguments to `recresid()`.


## strucchange 1.3-6

* Rnhanced documentation for new Rd parser.


## strucchange 1.3-5

* Added some further references to the vignette,
  and provide the associated .bib file in `~/inst/doc/`.

* Removed `\itemize` in Rd files for new R-devel.


## strucchange 1.3-4

* Fixed CITATION encoding.

* Removed Z.sty dependency in vignette.


## strucchange 1.3-3

* Enhanced references in the vignette, CITATION and 
  man pages.

* Fixed some outdated information in the vignette.


## strucchange 1.3-2

* Added new data set with bibliographic information about
  structural change publications.


## strucchange 1.3-1

* Renamed `SP500` to `SP2001` to avoid conflicts with `MASS`.


## strucchange 1.3-0

* Added NAMESPACE.

* Improved dependency declaration in DESCRIPTION.
