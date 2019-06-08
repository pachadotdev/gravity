
# gravity <img src="https://pachamaltese.github.io/gravity/hexicon.svg" width=150 align="right" alt="sticker"/>

[![Travis-CI Build
Status](https://travis-ci.org/pachamaltese/gravity.svg?branch=master)](https://travis-ci.org/pachamaltese/gravity)
[![cran
checks](https://cranchecks.info/badges/summary/gravity)](https://cran.r-project.org/web/checks/check_results_gravity.html)
[![Coverage
status](https://codecov.io/gh/pachamaltese/gravity/branch/master/graph/badge.svg)](https://codecov.io/github/pachamaltese/gravity?branch=master)
[![JOSS](http://joss.theoj.org/papers/10.21105/joss.01038/status.svg)](https://doi.org/10.21105/joss.01038)
[![Austrian Journal of
Statistics](https://zenodo.org/badge/doi/10.17713/ajs.v47i4.688.svg)](https://doi.org/10.17713/ajs.v47i4.688)

## About

Gravity models are used to explain bilateral flows related to the sizes
of bilateral partners, a measure of distance between them and other
influences on interaction costs.

The underlying idea is rather simple. The greater the masses of two
bodies and the smaller the distance between them, the stronger their
attraction. This concept is applied to several research topics such as
trade, migration or foreign direct investment. Even though their basic
idea is rather simple, gravity models can become complex regarding the
choice of models or estimation methods.

As especially for gravity beginners it is difficult to get an overview
of the different methods and implement them in R, the package provides a
wrapper of different standard estimation methods for gravity models. By
considering the descriptions and codes of the estimation methods, users
get a comprehensive and application-oriented overview of the different
methods, see which method may be suitable for a certain research
question or type of data and extend the code available to their research
projects.

The package provides estimation methods for log-log models and
multiplicative models.

Log-log estimation covers:

  - Ordinary Least Squares (OLS)
  - Fixed Effects
  - Double Demeaning (DDM)
  - Bonus vetus OLS with simple averages (BVU) and with GDP-weights
    (BVW)
  - Structural Iterated Least Squares (SILS)
  - Tetrads
  - Different Tobit models (Tobit, ET-Tobit, EK-Tobit)

Log-log models are partly complex to understand and program and thus a
comparison of them is not straightforward, wherefore the package aims at
easing an overview of the different methods combined with a direct
application.

Multiplicative estimation covers:

  - Poisson Pseudo Maximum Likelihood (PPML)
  - Gamma Pseudo Maximum Likelihood (GPML)
  - Negative Binomial Pseudo Maximum Likelihood (NBPML)
  - Nonlinear Least Squares (NLS)

Multiplicative models are covered with a log-link and different families
of distributions. These models are relatively easy to compute.

All functions named estimate gravity models, but they differ in whether
they estimate them in their additive or multiplicative form, their
requirements and assumptions with respect to data, their handling of
Multilateral Resistance terms as well as their possibilities concerning
the inclusion of unilateral independent variables. Therefore, they
normally lead to different estimation results. We refer the user to the
[Gravity Cookbook website](https://sites.google.com/site/hiegravity/)
for more information on gravity models in general.

Head and Mayer (2014) provide a comprehensive and accessible overview of
the theoretical and empirical development of the gravity literature as
well as the use of gravity models and the various estimation methods,
especially their merits and potential problems regarding applicability
as well as different datasets. In order to have a straightforward
application of all methods, the package comes with two example datasets,
one including and the other one excluding zero trade flows.

Examples utilizing these dataset are included in the description of the
methods. As the range of gravity models and their suitable estimation
methods is huge and specific to data and research question at hand, the
functions are kept simple such that researchers can utilize them as a
starting point of their research and get familiar with them.

On the [Gravity Cookbook
website](https://sites.google.com/site/hiegravity/), Keith Head and
Thierry Mayer provide Stata code for most methods. Where possible, the
functions are tested to lead the same results as the Stata code when
choosing the option of robust variance estimation. However, compared to
the Stata code available, the functions presented in this package
provide users with more flexibility regarding the type of estimation,
the number and type of independent variables as well as the possible
data.

Furthermore, additional functions are added. The functions of the
package include distance as an independent variable. If a panel model is
used time-invariant effects, such as geographical distance, may have to
be excluded; see the respective descriptions. Therefore, it is up to the
user to ensure that the functions can be applied to panel data.

Depending on the panel dataset and the variables - specifically the type
of fixed effects - included in the model, it may easily occur that the
model is not computable.

Depending on the specific model, the code of the respective function may
has to be changed in order to exclude the distance variable from the
estimation.

For a state-of-the-art exposition about cross-sectional data see Wölwer,
Breßlein, and Burgard (2018), and for a comprehensive overview of
gravity models for panel data see Egger and Pfaffermayr (2003),
Gómez-Herrera (2013) and Head, Mayer, and Ries (2010) as well as the
references therein (see also the references included in the descriptions
of the different functions).

At the very least, the user should take special care with respect to the
meaning of the estimated coefficients and variances as well as the
decision about which effects to include in the estimation. As, to our
knowledge at the moment, there is no explicit literature covering the
estimation of a gravity equation by Double Demeaning, Structural
Iterated Least Squares or Bonus vetus OLS using panel data, we do not
recommend to apply these methods in this case.

Contributions, extensions and error corrections are very welcome. Please
do not hesitate to contact us. We thank Martin Bresslein for valuable
comments and advice on early versions of the package.

## How to install

From R console:

``` r
# Install from CRAN
install.packages("gravity")

# Install from GitHub
devtools::install_github("pachamaltese/gravity")
```

## Community guidelines

If you want to contribute to the software, report issues or problems
with the software, please fork the repo and send us a Pull Request or
open an issue. We are happy to receive ideas and we would do our best to
coordinate efforts and improve this package without reinventing the
wheel.

If you seek support or have questions you can start a thread on the
issues section, or you can email us but we prefer open issues as
probably more users have the same questions as you.

## References

<div id="refs" class="references">

<div id="ref-Egger2003">

Egger, Peter, and Michael Pfaffermayr. 2003. “The Proper Panel
Econometric Specification of the Gravity Equation: A Three-Way Model
with Bilateral Interaction Effects.” *Empirical Economics* 28 (3):
571–80. <https://doi.org/10.1007/s001810200146>.

</div>

<div id="ref-Gomez-Herrera2013">

Gómez-Herrera, Estrella. 2013. “Comparing Alternative Methods to
Estimate Gravity Models of Bilateral Trade.” *Empirical Economics* 44
(3): 1087–1111. <https://doi.org/10.1007/s00181-012-0576-2>.

</div>

<div id="ref-Head2014">

Head, Keith, and Thierry Mayer. 2014. “Chapter 3 - Gravity Equations:
Workhorse,Toolkit, and Cookbook.” In *Handbook of International
Economics*, edited by Gita Gopinath, Elhanan Helpman, and Kenneth
Rogoff, 4:131–95. Handbook of International Economics. Elsevier.
<https://doi.org/10.1016/B978-0-444-54314-1.00003-3>.

</div>

<div id="ref-Head2010">

Head, Keith, Thierry Mayer, and John Ries. 2010. “The Erosion of
Colonial Trade Linkages After Independence.” *Journal of International
Economics* 81 (1): 1–14.
<https://doi.org/10.1016/j.jinteco.2010.01.002>.

</div>

<div id="ref-WoelwerBressleinBurgard2018">

Wölwer, Anna-Lena, Martin Breßlein, and Jan Pablo Burgard. 2018.
“Gravity Models in R.” *Austrian Journal of Statistics* 47 (4): 16–35.
<https://doi.org/10.17713/ajs.v47i4.688>.

</div>

</div>
