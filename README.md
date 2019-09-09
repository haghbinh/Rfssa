
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

This is a basic example which shows you how to use fssa algorithm to
analyze a univariate or multivariate functional time series:

## Univariate FSSA Example on Callcenter data

``` r
data("Callcenter")
require(fda)
require(Rfssa)
```

## Define functional objects

``` r
D <- matrix(sqrt(Callcenter$calls),nrow = 240)
N <- ncol(D)
time <- seq(ISOdate(1999,1,1), ISOdate(1999,12,31), by="day")
K <- nrow(D)
u <- seq(0,K,length.out =K)
d <- 22 #Optimal Number of basis elements
basis <- create.bspline.basis(c(min(u),max(u)),d)
Ysmooth <- smooth.basis(u,D,basis)
```

## Define functional time series

``` r
Y <- fts(Ysmooth$fd,time = time)
plot(Y,ylab = "Sqrt of Callcenter", xlab = "Intraday intervals")
```

## Univariate functional singular spectrum analysis

``` r
L <- 28
U <- fssa(Y,L)
plot(U,d=13)
```

![](man/figures/README-example4-1.png)<!-- -->

``` r
plot(U,d=9,type="lheats")
```

![](man/figures/README-example4-2.png)<!-- -->

``` r
plot(U,d=9,type="lcurves")
```

![](man/figures/README-example4-3.png)<!-- -->

``` r
plot(U,d=9,type="vectors")
```

![](man/figures/README-example4-4.png)<!-- -->

``` r
plot(U,d=10,type="periodogram")
```

![](man/figures/README-example4-5.png)<!-- -->

``` r
plot(U,d=10,type="paired")
```

![](man/figures/README-example4-6.png)<!-- -->

``` r
plot(U,d=20,type="wcor")
```

![](man/figures/README-example4-7.png)<!-- -->

``` r
gr <- list(1,2:3,4:5,6:7,8:20)
Q <- freconstruct(U, gr)
plot(Q[[1]],main="1st Component",ylab = " ", xlab = "Intraday intervals")
plot(Q[[2]],main="2nd Component",ylab = " ", xlab = "Intraday intervals")
plot(Q[[3]],main="3rd Component",ylab = " ", xlab = "Intraday intervals")
plot(Q[[4]],main="4th Component",ylab = " ", xlab = "Intraday intervals")
plot(Q[[5]],main="5th Component(Noise)",ylab = " ", xlab = "Intraday intervals")
```

### Other visiualisation types for object of class “fts”:

``` r
plot(Q[[1]], type="3Dsurface",
     main="1st Component",ylab = " ", xlab = "Intraday intervals")
```

![](man/figures/README-example5-1.png)<!-- -->

``` r
plot(Q[[2]][1:60], type="heatmap",
     main="2nd Component",ylab = " ", xlab = "Intraday intervals")
plot(Q[[3]][1:60], type = "3Dline",
     main="3rd Component",ylab = " ", xlab = "Intraday intervals")
```

![](man/figures/README-example5-2.png)<!-- --> \#\# Multivariate FSSA
Example on Bivariate Satelite Image Data

``` r
NDVI <- Jambi$NDVI
EVI <- Jambi$EVI
time <- Jambi$Date
```

## Kernel density estimation of pixel intensity

``` r
D0_NDVI <- matrix(NA,nrow = 512, ncol = 448)
D0_EVI <- matrix(NA,nrow =512, ncol = 448)
for(i in 1:448){
  D0_NDVI[,i] <- density(NDVI[,,i],from=0,to=1)$y
  D0_EVI[,i] <- density(EVI[,,i],from=0,to=1)$y
}
```

## Define functional objects

``` r
d <- 11
basis <- create.bspline.basis(c(0,1),d)
u <- seq(0,1,length.out = 512)
y_NDVI <- smooth.basis(u,as.matrix(D0_NDVI),basis)$fd
y_EVI <- smooth.basis(u,as.matrix(D0_EVI),basis)$fd
y <- list(y_NDVI,y_EVI)
```

## Define functional time series

``` r
Y <- fts(y,time=time)
plot(Y[1:100],ylab = c("NDVI","EVI"),main = "Probability Kernel Density")
```

## Multivariate functional singular spectrum analysis

``` r
L <- 45
U <- fssa(Y,L)
plot(U,d=10,type='values')
```

![](man/figures/README-example10-1.png)<!-- -->

``` r
plot(U,d=10,type='paired')
```

![](man/figures/README-example10-2.png)<!-- -->

``` r
plot(U,d=10,type='lheats', var = 1)
```

![](man/figures/README-example10-3.png)<!-- -->

``` r
plot(U,d=10,type='lcurves',var = 1)
```

![](man/figures/README-example10-4.png)<!-- -->

``` r
plot(U,d=10,type='lheats', var = 2)
```

![](man/figures/README-example10-5.png)<!-- -->

``` r
plot(U,d=10,type='lcurves',var = 2)
```

![](man/figures/README-example10-6.png)<!-- -->

``` r
plot(U,d=20,type='wcor')
```

![](man/figures/README-example10-7.png)<!-- -->

``` r

recon <- freconstruct(U = U, group = list(c(1),c(2,3),c(4)))
plot(recon[[1]],ylab = c("NDVI","EVI"), main = "Trend Components")
plot(recon[[2]],ylab = c("NDVI","EVI"), main = "Oscillation Components")
plot(recon[[3]],ylab = c("NDVI","EVI"), main = "Second Oscillation Components")
```
