#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: Yuchen Liu
Affiliation: Cavendish Astrophysics, University of Cambridge
Email: yl871@cam.ac.uk

Created in April 2023

Description: Posterior density plots based on the nested sampling chains
"""

import numpy as np
from anesthetic import read_chains, make_2d_axes
import matplotlib
import matplotlib.pyplot as plt
import os
from multiprocessing import Process

font = {'size':12}
matplotlib.rc('font', **font)

freq_boundaries = np.linspace(106, 196, 7, dtype=int)
freq_pair = [(freq_boundaries[i], freq_boundaries[i+1]) for i in range(6)]
os.system('mkdir -p posterior_plot')

def posterior_plot(pair):
    lower,upper = pair
    folder = 'ns_{}_{}MHz_estimate/gpr'.format(lower,upper)
    samples = read_chains(folder)
    # prior = samples.prior()
    names = ['p%s'%i for i in range(6)]
    fig, axes = make_2d_axes(names, figsize=(8, 8), facecolor='w',upper=False)
    samples.plot_2d(axes)
    plt.savefig('posterior_plot/posterior_gpr_model_{}_{}MHz.pdf'.format(lower,upper),bbox_inches='tight')
    
processes = []
for i in freq_pair:
    p = Process(target=posterior_plot, args=(i,))
    p.start()
    processes.append(p)
    
for p in processes:
    p.join()
