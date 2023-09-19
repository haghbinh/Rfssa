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
