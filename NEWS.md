Rfssa 1.0.0
===========

Updates
-------

-   `fssa` updated to allow for the decomposition stage of multivariate
    functional singular spectrum analysis (`mfssa`)

-   `freconstruct` updated to allow for the reconstruction stage of
    `mfssa`

-   `plot.fssa` was updated to allow the user to plot:

    1.  left singular functions in curves (`lcurves`)
    2.  left singular functions in heat map (`lheats`)
    3.  periodograms (`periodogram`)
    4.  right singular vectors (`vectors`)

-   `plot.fssa` was also updated to remove plotting options:

    1.  `meanvectors` and `meanpaired` removed in lieu of `vectors` and
        `paired` plot options
    2.  `efunctions` and `efunctions2` removed in lieu of `lcurves` and
        `lheats`

-   `wplot` was updated to allow the user to specify the `cuts`
    parameter to improve visualization of weighted correlation
    (w-correlation) matrix

New Additions
-------------

-   The functional time series (`fts`) class was added which extends
    functional data (`fd`) objects

-   `fts` arithmetic was added allowing the user to perform addition,
    subtraction, and multiplication of `fts` objects with other `fts`
    objects or `fts` objects with scalars

-   Indexing of `fts` objects is allowed

-   Correlation of `fts` objects can now be measured using `cor.fts`.

-   `plot.fts`, which uses plotly, was added to allow for visualization
    of `fts` objects

-   `launchApp` was added which allows the user to launch a shiny app
    for univariate functional singular spectrum analysis of multivariate
    functional singular spectrum analysis based on type specified. This
    app helps with visualization, and general exploration of the
    algorithm.

Minor improvements and bug fixes
--------------------------------

-   `Jambi` data set was added to create bivariate `mfssa` example

-   Plot ordering in `plot.fssa` was fixed for types `lcurves`,
    `lheats`, and `vectors`

-   More clear error checking was added to code
