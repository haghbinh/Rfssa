Rfssa 2.0.0
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

Minor improvements and bug fixes
--------------------------------

-   Fixed a bug in `line` type plots in `plot.fts` that prevented plotting.
