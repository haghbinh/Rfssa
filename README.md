

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Rfssa
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/Rfssa)](https://cran.r-project.org/package=Rfssa)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


The Rfssa package provides the collections of necessary functions to
implement functional singular spectrum analysis (FSSA)-based methods for 
analyzing univariate and multivariate functional time series (FTS). 
Univariate and multivariate FSSA are novel, non-parametric methods to perform decomposition and reconstruction of univariate and multivariate FTS 
respectively. In addition, the FSSA-based routines may be performed on FTS 
whose variables are observed over a one or two-dimensional domain. Finally, 
one may perform FSSA recurrent or vector forecasting of univariate or 
multivariate FTS observed over one-dimensional domains. Forecasting of FTS 
whose variables are observed over domains of dimension greater than one is 
under development.

# Introduction

The use of the package starts with the decomposition of functional time
series (fts) objects using the fssa routine. Then a suitable grouping of the
principal components is required for reconstruction (freconstruct) or 
forecasting (fforecast) which can be done heuristically by looking at the 
plots of the decomposition (plot). Once a suitable grouping is chosen, 
one may perform reconstruction where the sum of all the elements in the 
groups approximates the original FTS. One may also choose to perform 
forecasting after a grouping is chosen which returns future observations in 
the FTS.

# Updated Functionality

This version of the package leverages a new S4 object for FTS objects (fts). 
The new object may be specified using a provided basis and grid, a requested 
basis and grid, or a mixture of provided and requested elements. We note that 
the FTS object may be univariate or multivariate and variables may be observed 
over one or two-dimensional domains. Validity checking of the S4 object 
constructor inputs was also added to help guide the user. The plotting of FTS 
objects (fts.plot) was also updated to allow the user to plot FTS variables 
observed over two-dimensional domains. Next, the FSSA routine (fssa) was 
updated to perform faster by leveraging the RSpectra and RcppEigen R packages, 
and the Eigen C++ package. We achieved a roughly 20 times speed up for 
certain data examples. We updated the plotting of fssa objects (fssa.plot) to 
allow for plotting of left singular functions that correspond with FTS 
variables observed over a two-dimensional domain.
We updated FSSA reconstruction (freconstruct) to handle 
FTS whose variables are observed over one or two-dimensional domains. We also 
updated FTS arithmetic (such as FTS addition, FTS subtraction, etc.) to allow 
the user to perform scalar-FTS arithmetic on different variables of a 
multivariate FTS. In addition, we also now host the Callcenter, Jambi, and 
Montana datasets on GitHub to significantly decrease the size of the package. 
In order to load the data, one simply needs to use the load_github_data 
function. This same function can also be used to load data from any other 
public GitHub repository.

# New Functionality

The first piece of new functionality that has been added is that the user 
may now specify univariate or multivariate FTS comprised of variables observed 
over one or two-dimensional domains. In addition, forecasting of univariate 
and multivariate FTS observed over one-dimensional domains by FSSA/MFSSA 
recurrent forecasting and FSSA/MFSSA vector forecasting has also been added. 
We have also added in a new data set (Montana) which provides the data for a 
multivariate FTS observed over different dimensional domains.

# README Notes

The reader should note that we do not utilize FTS plotting options in
this README that are included in this update because of the large size
of the resulting files. The reader should refer to the ‘help’ of fssa to
see the same example but with the utilization of the new FTS plotting
options.

# Installation

You can install Rfssa from github with:

``` r
# install.packages("devtools")
devtools::install_github("haghbinh/Rfssa")
```

## Example

There are two basic examples which shows you how to use the 
basic fssa algorithm to analyze a univariate or multivariate functional time
series:

#### [FSSA (Univariate Visualization)](https://haghbinh.github.io/FSSA_report/uvisualization.html)
#### [FSSA (Univariate Decomposition)](https://haghbinh.github.io/FSSA_report/udecomposition.html)
#### [FSSA (Univariate Reconstruction)](https://haghbinh.github.io/FSSA_report/ureconstruction.html)
