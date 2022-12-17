# Automatic estimation of earthquake high-frequency strong-motion spectral decay

This code is the basis of a publication about automatic estimation of earthquake high-frequency strong-motion spectral decay (Sonnemann et al. 2019). It is provided as is without guarantees.

The goal of this study was to develop a flexible, efficient and automated algorithm for reliable and objective κ measurements from S wave spectra. The provided Matlab code automatically measures Kappa on the given accelerograms, using the method by Anderson & Hough 1984 of fitting a line to frequency vs. log(FAS), where FAS is the Fourier-domain amplitude spectrum of the S-wave signal from the acceleration record.

The P and S phase onset picker algorithm is not included, and these onset times are required input for each accelerometer record.

## How to run
- requires: Matlab, Signal Processing Toolbox (for `fir1()`)
    - Note: If the toolbox is unavailable, the filter function might be substituted by some alternative or omitted, however the change needs to be assessed.
- run `Kappa_demo_sa.m`
    - loads example dataset `MyData.mat`
    - runs main function of `K_automeasure_dyn2d_sa.m`
- `K_automeasure_dyn2d_sa.m` does all the work, returns output and can plot
    - output values: Kappa of each component, various diagnostics
    - plots if input `FIG='on'`
- details: please read the function descriptions at the top

## References
- Sonnemann, T., Halldorsson, B., & Jonsson, S. (2019). Automatic estimation of earthquake high-frequency strong-motion spectral decay in south Iceland. Soil Dynamics and Earthquake Engineering, 125, 105676. [DOI](https://doi.org/10.1016/j.soildyn.2019.05.015)

