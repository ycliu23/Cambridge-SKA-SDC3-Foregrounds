# Camrbidge SKA Science Data Challenge - Foregrounds
## OSKAR Simulation
The SKA End-to-End simulation pipeline used in the SKA SDC3 can be found at
https://github.com/ycliu23/SKA_Power_Spectrum_and_EoR_Window, 
which is a branch of https://github.com/oharao/SKA_Power_Spectrum_and_EoR_Window with a modified sky model (inner sky + outer sky). We use the GLEAM and LOBES catalogue provided to simulate a 4-hour tracking observation of point sources that consists of 1440 time steps.
(point source calibration)
## Imaging
The imaging process utilizes WSCLEAN (Offringa et al., 2014) to grid and inverse Fourier transform the visibilities from OSKAR to images by properly accounting for the sky curvature (ie. w-terms).
## Foreground Removal
## MCMC Sampling
## Power Spectrum Analysis
