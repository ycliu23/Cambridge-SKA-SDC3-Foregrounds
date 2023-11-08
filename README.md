# Camrbidge SKA Science Data Challenge - Foregrounds
## OSKAR Simulation
The SKA End-to-End simulation pipeline used in the SKA SDC3 can be found at https://github.com/ycliu23/SKA_Power_Spectrum_and_EoR_Window.git, which is a branch of https://github.com/oharao/SKA_Power_Spectrum_and_EoR_Window with a modified sky model (inner sky + outer sky).
(point source calibration)
## Imaging
The imaging process utilizes WSCLEAN (Offringa et al., 2014) to grid and inverse Fourier transform the visibilities from OSKAR to images by properly accounting for the sky curvature (ie. w-terms).
## Foreground Removal
## MCMC Sampler
## Power Spectrum Analysis
