<img src="https://zeileis.codeberg.page/strucchange/strucchange-icon.png" align="right" alt="strucchange icon" width="100" />

# Testing, Monitoring, and Dating Structural Changes


## Overview

The R package [strucchange](https://zeileis.codeberg.page/strucchange/) provides
a comprehensive toolbox for testing, monitoring, and dating structural changes in
linear regression models. Many of the methods have also been generalized to any
parametric model estimated by least squares, maximum likelihood, and other M-type
estimators. In short, these methods are concerned with answering the following questions.

* _Testing:_ Are the parameters of a model stable throughout the sample period or
  is there evidence that they changed over time?
* _Monitoring:_ If a model with stable parameters could be established, do the
  parameters remain stable as new observations come in?
* _Dating:_ If there is evidence for changes in the parameters, when and how did the
  parameters change?

Various families of tests are implemented, including the 
generalized fluctuation test framework as well as
the $F$ test (or Chow test) framework. This includes methods to
fit, plot and test fluctuation processes (e.g., CUSUM, MOSUM,
recursive/moving estimates) and $F$ statistics, respectively.


## Citations

Zeileis A, Leisch F, Hornik K, Kleiber C (2002).
  "strucchange: An R Package for Testing for Structural Change in Linear Regression Models."
  _Journal of Statistical Software_, **7**(2), 1-38.
  [doi:10.18637/jss.v007.i02](https://doi.org/10.18637/jss.v007.i02)

Zeileis A, Kleiber C, Krämer W, Hornik K (2003).
  "Testing and Dating of Structural Changes in Practice."
  _Computational Statistics & Data Analysis_, **44**(1-2), 109-123.
  [doi:10.1016/S0167-9473(03)00030-6](https://doi.org/10.1016/S0167-9473%2803%2900030-6)

Zeileis A (2006).
  "Implementing a Class of Structural Change Tests: An Econometric Computing Approach."
  _Computational Statistics & Data Analysis_, **50**(11), 2987-3008.
  [doi:10.1016/j.csda.2005.07.001](https://doi.org/10.1016/j.csda.2005.07.001)



## Installation

The stable version of `strucchange` is available from
[CRAN](https://CRAN.R-project.org/package=strucchange):

``` r
install.packages("strucchange")
```

The latest development version can be installed from
[R-universe](https://zeileis.R-universe.dev/strucchange):

``` r
install.packages("strucchange", repos = "https://zeileis.R-universe.dev")
```


## License

The package is available under the
[General Public License version 3](https://www.gnu.org/licenses/gpl-3.0.html)
or [version 2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
