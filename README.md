

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Rfssa
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/Rfssa)](https://cran.r-project.org/package=Rfssa)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


The Rfssa package provides the collection of necessary functions to
implement functional singular spectrum analysis (FSSA)-based methods for
analyzing univariate and multivariate functional time series (FTS).
Univariate and multivariate FSSA are novel, non-parametric methods used to perform decomposition and reconstruction of univariate and multivariate FTS respectively. In addition, the FSSA-based routines may be performed on FTS whose variables are observed over a one or two-dimensional domain. Finally, one may perform FSSA recurrent or FSSA vector forecasting of univariate or multivariate funts observed over one-dimensional domains. Forecasting of funts whose variables are observed over domains of dimension greater than one is under development.

# Summary

The use of the package begins by defining an `funts' object by providing the constructor with the raw data, basis specifications, and grid specifications. We note that the FTS object may be univariate or multivariate and variables may be observed
over one (curves) or two-dimensional (images) domains. Validity checking of the S4 object constructor inputs is included to help guide the user. The user may leverage the plot.funts method to visualize the funts object. A variety of plotting options are available for variables observed over a one-dimensional domain and a visuanimation is offered for variables observed over a two-dimensional domain. Next, the user provides the funts object and a chosen lag parameter to the FSSA routine (fssa) to obtain the decomposition. We note that the decomposition function leverages the RSpectra and RcppEigen R packages, and the Eigen C++ package to speed up the routine. The plot.fssa method may be used to visualize the results of the decomposition and to choose an appropriate grouping of the eigentriples for reconstruction (freconstruct) or forecasting (fforecast). The freconstruct routine can be used to reconstruct a list of funts objects specified by the grouping while the fforecast function returns a list of funts objects that contain predictions of the signals specified by the grouping. The user may also calculate the bootstrapped prediction interval for forecasts using the fpredinterval function. We note that when forecasting is performed, usually the user specifies one group that captures the assumed deterministic, extracted signal that is found within the FTS and all other modes of variation are excluded. We also note that currently, forecasting only supports FTS whose variables are observed over a one-dimensional domain with two-dimensional domain forecasting to be added in the future.

Other functionalities offered by the package include:
<ul>
  <li> funts arithmetic - Allows the user to perform FTS-FTS arithmetic and FTS-scalar arithmetic (such as addition, subtraction, etc.).</li>
  <li> load_github_data - Allows the user to load any .RData file hosted on GitHub including the Callcenter, Jambi,
and Montana datasets.</li>
  <li> launchApp - Launches the built-in R Shiny app that can be used to interactively explore the FSSA-based routines on various datasets.</li>
</ul>

# Updates

- The name `fts` has been modified to `funts` to avoid any clashes with the  package. Furthermore, the class of `funts` has bee transitioned from S4 to S3 to ensure better compatibility and consistency within the package.
These changes are aimed at preventing any conflicts when using `Rfssa` in conjunction with other packages like `rainbow`, enhancing the user experience.

-  All the methods for `funts` have re-implemented and introduced new generic methods such as `length()`, `print()`, and `plot()` to provide a more comprehensive and user-friendly interface. 

- The `plot()` method for `funts` class objects (formerly `fts`) has been renamed to `plotly_funts()`. This new name more accurately reflects the type of plots it generates, which are based on `plotly` graphics.

- An S3 class named `fforecast` is added to encapsulate the output of the `fforecast()` function. This class is designed to provide a more organized and intuitive structure for handling forecasted functional time series (FTS) data.

- Three convenient functions, namely `loadJambiData()`, `loadCallcenterData()`, and `loadMontanaData()` are added. These functions have been designed to simplify the process of acquiring the raw dataset from the web and loading it into the global environment.
    
- In the latest version of the package, two new parameters, start and end, have been introduced in the `funts` function to capture the duration of the time series. These parameters provide flexibility for users to specify time information in a more structured and standardized manner. Users can now set `start` and `end` using various time and date classes such as `Date`, `POSIXct`, or `POSIXt`, allowing for better representation of time.


# README Notes

The reader should note that we do not utilize FTS plotting options in
this README because of the large size
of the resulting files. The reader should refer to the examples offered at the end of this README 
to see examples of how to apply the methodologies to real data.

# Installation

You can install Rfssa from github with:

``` r
# install.packages("devtools")
devtools::install_github("haghbinh/Rfssa")
```

# Examples

The following links provide examples of how to run FSSA-related methods on real data:


### [FSSA Visualization](https://jtrinka.github.io/visualization.html)
### [FSSA Decomposition](https://jtrinka.github.io/decomposition.html)
### [FSSA Reconstruction](https://jtrinka.github.io/reconstruction.html)
### [FSSA Forecasting](https://jtrinka.github.io/forecasting.html)

