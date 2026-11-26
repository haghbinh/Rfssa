 Rfssa 3.2.1.9001 (development version, branch AESSA)
===========

New features (experimental)
-------------

- **Autoencoder SSA (AESSA family).** New functions `aessa()`, `vaessa()`, `aeforecast()`,
  `aessa_select()`, `aessa_arch()` and a `reconstruct()` method replace the rank-r SVD step
  of SSA by an autoencoder of the lagged records (nonlinear SSA).
  - `aessa()`: deterministic autoencoder; `vaessa()`: variational version with Gaussian
    decoder, isotropic or per-series (block) noise, a diagonal-Gaussian amortized encoder and
    a composite ELBO; the nonlinear counterpart of `prssa()`.
  - Default architecture: decoder = linear part + tanh residual network whose output layer
    starts at zero, with the linear part initialized at Basic SSA (`aessa()`) or PrSSA
    (`vaessa()`), so the untrained model is exactly SSA/PrSSA. Restarts, early stopping on a
    validation block of windows separated by a gap, refit on all windows, seeds,
    float64/float32.
  - Inputs: numeric vector/`ts`, numeric matrix (multivariate, vertical stacking as in
    PrMSSA), `funts` (Gram-whitened basis coefficients as in PrFSSA) and lists of these.
  - Forecasting (`aeforecast()`): per-horizon decoder inversion (window conditioning,
    Levenberg-Marquardt with multi-start), recursive inversion, or a VAR(1) on the latent
    codes; VAESSA predictive intervals by a Laplace approximation (exact for a linear
    decoder, equal to the PrSSA window law) or Monte Carlo.
  - Latent-dimension selection by rolling-origin validation (`aessa_select()`).
  - Methods `print()`, `summary()`, `plot()` (loss, latent codes, components), `predict()`.
- The training backend is the `torch` package, which is suggested, not required; the AESSA
  functions stop with an informative message when `torch` (or LibTorch) is not installed.
  All other functions are unchanged.

 Rfssa 3.2.1.9000 (development version)
===========

New features
-------------

- **Probabilistic SSA (PrSSA family).** New functions `prssa()`, `prreconstruct()`,
  `prforecast()` and `prrank()` fit a Gaussian latent-factor model to the lagged records
  of a series by composite likelihood. They cover PrSSA (univariate), PrMSSA (several series,
  vertical stacking), PrFSSA (`funts` objects, exact L2 geometry through the Gram matrix of
  the basis) and hybrid PrHSSA (lists mixing numeric series and `funts` objects).
  - Noise models: isotropic (closed form; with `mean = "zero"` exactly Basic SSA/FSSA/VMSSA),
    one level per series or variable (ECME, i.e. weighted MSSA with self-consistent
    weights), and heteroscedastic functional noise `log sigma^2(s) = b(s)'eta` with a
    P-spline penalty, for `funts` variables or for a numeric matrix of curves on a grid
    (rows = times, columns = grid points; argument `grid`).
  - Smoothing by PPCA of half-smoothed records (`alpha > 0`, regularized HSSA in the
    noise-free limit), with a second-stage raw-scale noise level for prediction.
  - Reconstruction by posterior means (shrunken SSA) or projections.
  - Forecasting by the window-conditional law (h <= L - r; per-horizon marginal laws by
    default, each horizon conditioned on its own window, or with `joint = TRUE` the joint
    law of all window horizons and their covariance `joint_cov`), a VAR(L-1) recursion with
    variance propagation beyond, or vector mode (posterior-mean seed, then the noise-free
    V-forecast operator of the fitted loadings, so the shrinkage is applied once); pointwise
    and simultaneous prediction bands (functional bands for `funts`).
  - Rank selection by nominal CL-BIC (default; the most accurate in the simulation study of
    the PrSSA paper), CL-BIC with a sandwich (HAC) penalty, or CL-AIC.
  - Methods `print()`, `summary()`, `logLik()` and `plot()` for `prssa` fits.
- Added testthat tests of the exact reductions to Rssa (Basic SSA, weighted MSSA, `lrr`,
  `rforecast`, `vforecast`) and to `fssa()` with an orthonormal basis.
- Existing functions are unchanged.



 Rfssa 3.2.1
===========

Changes
-------------

- **CRAN policy on Internet access:** `loadCallcenterData()`, `loadJambiData()`,
  `loadMontanaData()`, `loadAustinData()` and `loadUtqiagvikData()` now fail gracefully.
  If a dataset cannot be downloaded from GitHub (no Internet connection, timeout, or a
  changed resource), they return `NULL` invisibly with an informative message instead of
  an error. Their examples run inside `\donttest{}` and check for `NULL`. The `funts()`
  example now also includes an offline example with simulated data, and the Shiny apps
  (`launchApp()`) show a message instead of stopping when the example data cannot be
  downloaded. (Version 3.2.0 was archived on CRAN on 2026-04-28 for this reason.)

- `plot.funts()` now restores the user's graphical parameters (`par()`) on exit.



 Rfssa 3.2.0
===========

Changes
-------------

- **Improved Rd Cross-References:** 
  Fixed `\link{}` targets in documentation by adding explicit package anchors (e.g., `\link[fda]{fd}`), in compliance with updated CRAN HTML reference manual policies.

- 🛠️ **Build Compatibility:**
  Ensured compatibility with current R build toolchains by verifying successful compilation on Windows using Rtools.



Rfssa 3.1.0
===========

New Additions
-------------
 The new following functions are replaced to GitHub load data:
 
-  `loadCallcenterData()`
-  `loadJambiData()`
-  `loadMontanaData()`
-  `loadAustinData()`
-  `loadUtqiagvikData()`

In these updated functions, upon downloading the data files from GitHub into a temporary directory (not the global environment), the target objects are now returned within the function. This modification allows users to save the data into an arbitrary variable of their choice.


Rfssa 3.0.2
===========

Minor improvements and bug fixes
--------------------------------
-   Fixed a bug in title of the `3Dline` and `3DSurface` type plots in `plotly_funts` that prevented plotting.





Rfssa 3.0.0
===========

Updates
-------

- The name `fts` has been modified to `funts` to avoid any clashes with the  package. Furthermore, the class of `funts` has bee transitioned from S4 to S3 to ensure better compatibility and consistency within the package.
These changes are aimed at preventing any conflicts when using `Rfssa` in conjunction with other packages like `rainbow`, enhancing the user experience.

-  All the methods for `funts` have re-implemented and introduced new generic methods such as `length()`, `print()`, and `plot()` to provide a more comprehensive and user-friendly interface. 

- The `plot()` method for `funts` class objects (formerly `fts`) has been renamed to `plotly_funts()`. This new name more accurately reflects the type of plots it generates, which are based on `plotly` graphics.

New Additions
-------------

- An S3 class named `fforecast` is added to encapsulate the output of the `fforecast()` function. This class is designed to provide a more organized and intuitive structure for handling forecasted functional time series (FTS) data.

- Three convenient functions, namely `loadJambiData()`, `loadCallcenterData()`, and `loadMontanaData()` are added. These functions have been designed to simplify the process of acquiring the raw dataset from the web and loading it into the global environment.
    
Minor improvements and bug fixes
--------------------------------

- In the latest version of the package, two new parameters, start and end, have been introduced in the `funts` function to capture the duration of the time series. These parameters provide flexibility for users to specify time information in a more structured and standardized manner. Users can now set `start` and `end` using various time and date classes such as `Date`, `POSIXct`, or `POSIXt`, allowing for better representation of time.




Rfssa 2.0.1
===========

Updates
-------
-   `fts` updated to be a constructor of a new S4 object that is used to create
    objects of class `fts`. Note that the user may specify their own basis 
    and grid or may specify that they want the constructor to build such 
    attributes automatically. In addition, the `fts` may be comprised of 
    variables observed over one or two-dimensional domains. This constructor 
    also has custom validity checks built in to help the user construct `fts` 
    objects.
    
-   `fts.plot` was updated to allow for plotting of `fts` variables observed 
    over two-dimensional domains.

-   `fssa` updated so that univariate and multivariate `fts` objects may be 
    decomposed where each variable may be observed over a one or 
    two-dimensional domain. In addition, the speed of the decomposition 
    process was increased by using the Rspectra and RcppEigen R packages, 
    and the Eigen C++ package.
    
-   `plot.fssa` was updated to allow the user to plot the left singular 
    functions for `fts` whose variables have domains over two-dimensions. All 
    other plotting options were updated to handle the two-dimensional domain 
    functionality as well.

-   `freconstruct` updated to allow for the reconstruction stage of
    univariate and multivariate `fts` objects whose variables might be 
    might observed over one or two-dimensional domains.

-   Arithmetic operations such as `fts` addition and `fts` subtraction have 
    been updated to allow for numeric vector-`fts` arithmetic.


New Additions
-------------

-   `fforecast` was added to allow for nonparametric forecasting of `fts` 
    objects via `fssa` recurrent or `fssa` vector forecasting. The 
    approach begins with objects of class `fssa`. The `fts` may be 
    univariate or multivariate however, the variables must be observed over 
    a one-dimensional domain. Different dimensional domains forecasting is 
    under development.
    
-   `Montana` was added to provide an example of a multivariate `fts` 
    whose variables are observed over a one-dimensional and two-dimensional 
    domain.
    
-   `load_github_data` was added to allow the user to load the `Callcenter`, 
    `Jambi`, and `Montana` datasets from GitHub which significantly reduced 
    the size of the package. In addition, this same function may be used to 
    load data from any other public GitHub repository.

Minor improvements and bug fixes
--------------------------------

-   Fixed a bug in `line` type plots in `plot.fts` that prevented plotting.

