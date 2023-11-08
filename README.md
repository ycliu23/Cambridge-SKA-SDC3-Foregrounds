# Camrbidge SKA Science Data Challenge - Foregrounds
## OSKAR Simulation
The SKA End-to-End simulation pipeline used in the SKA SDC3 can be found at
https://github.com/ycliu23/SKA_Power_Spectrum_and_EoR_Window, 
which is a branch of https://github.com/oharao/SKA_Power_Spectrum_and_EoR_Window with a modified sky model (inner sky + outer sky). We use the GLEAM and LOBES catalogue provided by the SKAO to simulate a 4-hour track observation of point sources that consists of 1440 time steps, each integrating over 10 seconds. The observation also covers the same frequency range as in the SDC3 from 106 MHz to 196 MHz, with a interval of 0.1 MHz.
## Imaging
The imaging process utilizes WSCLEAN (Offringa et al., 2014) to grid and inverse Fourier transform the visibilities from OSKAR to images by properly accounting for the sky curvature (ie. w-terms).
## Foreground Removal
## MCMC Sampling
## Power Spectrum Analysis
