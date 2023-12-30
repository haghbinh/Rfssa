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

