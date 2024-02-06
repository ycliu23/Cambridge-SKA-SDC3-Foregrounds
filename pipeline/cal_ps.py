#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: Yuchen Liu
Affiliation: Cavendish Astrophysics, University of Cambridge
Email: yl871@cam.ac.uk

Created in April 2023

Description: Power spectral analysis of foreground-subtracted residual visibilities
"""

from ps_eor import datacube, pspec, psutil
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import astropy.constants as const
from scipy import stats
import os
from multiprocessing import Process
import sys
import re

font = {'size':12}
matplotlib.rc('font', **font)

inpath = 'gpr_processed_cube/'
outpath = 'figures/'
os.system('mkdir -p {}'.format(outpath))

FREQ_TAPER = 'nuttall'
ERROR_P1D = 1
PREFIX = 'gpr_res'
MASK = True
THRESHOLD = 9

K_PERP = np.linspace(5e-2,5e-1,10)
K_PAR = np.linspace(5e-2,5e-1,10)

def ratio(array):
    ratios = np.zeros(len(array)-1)
    for i in range(len(ratios)):
        ratios[i] = array[i] / array[i+1]
    return ratios

def error_prop(array):
    return np.sqrt(np.sum(array**2))

res_cubes = []
for i in os.listdir(inpath):
    if i.endswith('.h5'):
        res_cubes.append(i)
res_cubes = np.sort(res_cubes)

def cal_ps(file):
    freq_lower, freq_upper, _ = np.array(re.findall(r'\d+',file)).astype(int)
    kind = file.split('_')[1]
                                                                    
    data_cube = datacube.CartDataCube.load(inpath+file)
        
    eor_bin_list = pspec.EorBinList(data_cube.freqs)
    eor_bin_list.add_freq(1,freq_lower,freq_upper)
    eor = eor_bin_list.get(1)
    el = psutil.k_to_l(K_PERP,eor.z)
    ps_conf = pspec.PowerSpectraConfig(el,FREQ_TAPER)
    ps_conf.umax = data_cube.ru.max()
    ps_conf.umin = data_cube.ru.min()
    pb = datacube.SkaLowPrimaryBeam()
    ps_gen = pspec.PowerSpectraCart(eor, ps_conf, pb)

    k_par_centre = K_PAR
    dk_par_centre = np.diff(k_par_centre)
    k_par_bins = np.concatenate([[k_par_centre[0]-dk_par_centre[0]/2],k_par_centre[:-1]+dk_par_centre/2,[k_par_centre[-1]+dk_par_centre[-1]/2]])
    
    ps2d = ps_gen.get_ps2d(data_cube)
    k_par = ps2d.k_par
    k_per = ps2d.k_per
    ps_data = ps2d.data
    ps_err = ps2d.err
    binned_data = np.zeros((len(k_par_centre),len(k_per)))
    binned_error = np.zeros((len(k_par_centre),len(k_per)))
    
    for i in range(len(k_per)):
        binned_data[:,i] = stats.binned_statistic(k_par,ps_data[:,i],bins=k_par_bins).statistic
        binned_error[:,i] = stats.binned_statistic(k_par,ps_err[:,i],bins=k_par_bins,statistic=error_prop).statistic
    assert binned_data.shape == binned_error.shape
    print(binned_data.shape)
    
    if MASK:
        k_per,k_par_centre,binned_data,binned_error
        for icol in range(binned_data.shape[1]):
            col = binned_data[:,icol]
            col_ratios = ratio(col)
            i_mask = np.where(col_ratios > THRESHOLD)[0][-1]
            if (icol < 5) & (i_mask > 2):
                i_mask = 2
            binned_data[:,icol][:i_mask+1] = np.nan
            binned_error[:,icol][:i_mask+1] = np.nan

    plt.figure()
    plt.pcolormesh(k_per,k_par_centre,binned_data,norm='log',cmap='jet')
    plt.xlim(k_per.min(),k_per.max())
    plt.ylim(k_par_centre.min(),k_par_centre.max())
    plt.xscale('log')
    plt.yscale('log')
    plt.colorbar()
    plt.savefig(outpath+'ps2d_res_{}_{}MHz'.format(freq_lower,freq_upper)+'.png',bbox_inches='tight')
    
    plt.figure()
    plt.pcolormesh(k_per,k_par_centre,binned_error,norm='log',cmap='jet')
    plt.xlim(k_per.min(),k_per.max())
    plt.ylim(k_par_centre.min(),k_par_centre.max())
    plt.xscale('log')
    plt.yscale('log')
    plt.colorbar()
    plt.savefig(outpath+'ps2d_err_res_{}_{}MHz'.format(freq_lower,freq_upper)+'.png',bbox_inches='tight')
    
    ### Submission
    np.savetxt(outpath+'Cantabrigians_{:.1f}MHz-{:.1f}MHz.data'.format(freq_lower,freq_upper),binned_data)
    np.savetxt(outpath+'Cantabrigians_{:.1f}MHz-{:.1f}MHz_errors.data'.format(freq_lower,freq_upper),binned_error)
    
processes = []
for i in res_cubes:
    p = Process(target=cal_ps, args=(i,))
    p.start()
    processes.append(p)
    
for p in processes:
    p.join()
