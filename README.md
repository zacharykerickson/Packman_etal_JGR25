# Packman_etal_JGR25
Code to conduct analyses and make figures found in Packman et al. (submitted to JGR-Oceans)

## Authorship
The code contained herein was first written by Sarah Packman (Harvard) during a NOAA Lapenta Internship at NOAA Pacific Marine Environmental Laboratory (PMEL), and afterwards modified by Zachary Erickson (PMEL).

## Python and Jupyter files:

- trend_data.py - inputs RFROM .nc files and outputs the seasonal trend coefficients at distinct pressure levels into the file trend.nc  
- global_analyses_variability.ipynb - inputs RFROM .nc files and trend.nc to detrend data and produce metrics for temperature anomaly standard deviations, which is output as stdev_T.nc  
- correlation_analysis_final.ipynb - inputs RFROM .nc files and trend.nc to detrend data and perform a lagged correlation analysis, which is output as lagged_correlation_analysis_150m_final.nc  
- Figures_final.ipynb - inputs RFROM .nc files, stdev_T.nc, and lagged_correlation_analysis_150m_final.nc to reproduce the 6 figures in the paper (and 2 supplemental figures), which are saved as .png files

## netCDF files

- stdev_T.nc - output from global_analyses_variability.ipynb. Contains the data used to make Figure 2 in the paper.
  
Note: trend.nc (160 MB) and lagged_correlation_analysis_150m_final.nc (243 MB) are not saved here due to their size, but can be reconstructed using RFROM .nc files, found at https://data.pmel.noaa.gov/pmel/erddap/files/argo_rfromv21_temp/

## Figures

- Figure_timeseries_final.png - Figure 1 of paper  
- Figure_stdevT_final.png - Figure 2 of paper  
- Figure_anomalies_final.png - Figure 3 of paper  
- Figure_lagged_correlations_final.png - Figure 4 of paper  
- Figure_correlation_maps_final.png - Figure 5 of paper
- Figure_map_types_final.png - Figure 6 of paper
- Figure_lagged_correlation_supplement_final.png - Figure S1 of paper
- Figure_lagged_correlation_omit_supplement_final.png - Figure S2 of paper
