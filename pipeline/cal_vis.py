#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: Yuchen Liu
Affiliation: Cavendish Astrophysics, University of Cambridge
Email: yl871@cam.ac.uk

Created in April 2023

Description: Fourier transform image cube to gridded visibilities followed by PSF deconvolution
"""

from ps_eor import datacube, pspec
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import astropy.constants as const
import os
from multiprocessing import Process

font = {'size':12}
matplotlib.rc('font', **font)

path = ''

FOV = 4
SPATIAL_TAPER = 'tukey'
FREQ_TAPER = 'nuttall'
ERROR_P1D = 1
WEIGHT = ['msn'] # ['msn','msw']

def cal_ps(sub_folder,weighting):
    inpath = path + '/' + sub_folder + '/'
    freq_lower = int(sub_folder.split('_')[0])
    freq_upper = int(sub_folder.split('_')[1].rstrip('MHz'))
    freq_central = np.mean([freq_lower,freq_upper])
    b_lambda_min = 30
    b_lambda_max = 73.5e3/(const.c.value/(freq_central*1e6))
    for i in os.listdir(inpath):
        if i.endswith('image.fits'):
            if weighting == i.split('_')[2]:
                image = i
        elif i.endswith('psf.fits'):
            psf = i
            
    if SPATIAL_TAPER == None:
        spatial_window = None
    else:
        spatial_window = datacube.WindowFunction(SPATIAL_TAPER)
        
    data_cube = datacube.CartDataCube.load_from_fits_image_and_psf([inpath+image],[inpath+psf],
                                                                b_lambda_min, b_lambda_max, np.radians(FOV),compat_wscnormf='new',window_function=spatial_window,
                                                                int_time=10,total_time=1000,min_weight_ratio=0.02)
                                                                    
    data_cube.save(inpath+'processed_cube_{lower}_{upper}MHz_{weighting}.h5'.format(lower=freq_lower,upper=freq_upper,weighting=weighting))
    

sub_bands = os.listdir(path)
combination = [[i, j] for i in sub_bands for j in WEIGHT]
processes = []
for i in combination:
    p = Process(target=cal_ps, args=i)
    p.start()
    processes.append(p)
    
for p in processes:
    p.join()
