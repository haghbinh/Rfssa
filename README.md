
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Rfssa

The Rfssa package provides the collections of necessary functions to
implement Functional Singular Spectrum Analysis (FSSA) for analysing
Functional Time Series (FTS). FSSA is a novel non-parametric method to
perform decomposition and reconstruction of FTS.

# Introduction

The use of the package starts with the decomposition of functional time
series (fts) objects using fssa. Then a suitable grouping of the
principal compomnents is required for reconstruction which can be done
heuristically by looking at the plots of the decomposition (plot).
Alternatively, one can examine the weighted correlation (w-correlation)
matrix (fwcor). The final step is the reconstruction of the principal
components into additive fts objects whose sum approximates the original
univariate or multivariate functional time series (freconstruct).

# Updated Functionality

This version of the Rfssa package includes updates to existing functions
including fssa, plot, wplot, and freconstruct. Multivariate functional
singular spectrum analysis (mfssa) was added to the package in fssa to
allow the user to perform embedding and decomposition of a multivariate
FTS. The reconstruction stage in freconstruct was also updated to allow
for reconstruction (including Hankelization) of multivariate FTS objects
using multivariate FSSA objects that come from mfssa. Plotting options
for FSSA objects in plot were also updated so that the user can now plot
left singular functions, right singular vectors, left singular function
heat diagrams, and periodograms. FSSA plotting options also allow the
user to specify which particular components they want to plot. For
example, a user can specify that they want to see a paired-plot of only
the third and fourth component. The ‘meanvectors’ and ‘meanpaired’
options were removed as these are satisfied with ‘paired’ and ‘vectors’
options. The ‘efunctions’ and ‘efunctions2’ options were also removed in
lieu of the addition of the left singular function heat map option. The
user can also specify the ‘cuts’ parameter in wplot to make
visualization of the w-correlation matrix easier.

# New Functionality

This version of the Rfssa package also includes new functions for
converting functional data (FD) objects to FTS objects, arithmetic,
indexing, correlation, and plotting of FTS data. The user is able to
convert an FD object to an FTS object using fts. The user can also
perform addition, subtraction, and multiplication of FTS objects with
other FTS objects or FTS objects with scalars by using ‘+’, ‘-’, and
’\*’ respectively. The package also allows for indexing of FTS
objects by using ‘\[’. The user can also measure the unweighted
correlation between FTS objects by using cor.fts. The plotting of FTS
objects can be performed using plot which uses the plotly package for
visualization.

The package update also includes a new shiny app (launchApp) that can be
used for demonstrations of univariate or multivariate FSSA depending on
the type that is specified. The app allows the user to explore FSSA with
simulated data, data that is provided on the server, or data that the
user provides. It allows the user to change parameters as they please,
gives visual results of the methods, and also allows the user to compare
FSSA results to other spectrum analysis methods such as multivariate
singular spectrum analysis. The tool is easy to use and can act as a
nice starting point for a user that wishes to perform FSSA as a part of
their data analysis.

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

There are two basic examples which shows you how to use fssa algorithm
to analyze a univariate or multivariate functional time
series:

### [FSSA (Part I: Decomposition Step)]<https://haghbinh.github.io/FSSA_report/Jambi2.html>

### [FSSA (Part II: Reconstruction Step)]<https://haghbinh.github.io/FSSA_report/Jambi3.html>

### [MSSA Comparing with FSSA]<https://haghbinh.github.io/FSSA_report/Jambi4.html>
