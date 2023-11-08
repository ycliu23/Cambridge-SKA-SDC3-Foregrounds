from ps_eor import datacube, pspec, psutil
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import astropy.constants as const
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
    el = 2 * np.pi * (np.arange(data_cube.ru.min(), data_cube.ru.max(), 10))
    ps_conf = pspec.PowerSpectraConfig(el,FREQ_TAPER)
    ps_conf.umax = data_cube.ru.max()
    ps_conf.umin = data_cube.ru.min()
    pb = datacube.SkaLowPrimaryBeam()
    ps_gen = pspec.PowerSpectraCart(eor, ps_conf, pb)
    print(ps_gen.config.umax)
    print(ps_gen.config.umin)
    
    ps2d = ps_gen.get_ps2d(data_cube)
    k_par = ps2d.k_par
    k_per = ps2d.k_per
    ps_data = ps2d.data
    ps_err = ps2d.err
    
    ix, iy = np.where(~np.isnan(ps_data))
    iy = np.unique(iy)
    ps_data = ps_data[:,iy]
    ps_err = ps_err[:,iy]
    k_per = k_per[iy]
    # print(ps_data.max())
    # print(ps_data.min())
    
    freq_inteval = str(freq_lower)+'_'+str(freq_upper)+'MHz'
    
    ps2d.save_to_txt(outpath+'ps2d_'+freq_inteval+'.txt')
    plt.figure()
    plt.pcolormesh(k_per,k_par,ps_data,norm='log',cmap='jet')
    plt.xlim(k_per.min(),k_per.max())
    plt.ylim(k_par.min(),k_par.max())
    plt.xlabel('$k_\perp\ [h\ \mathrm{cMpc^{-1}}$]')
    plt.ylabel('$k_\parallel\ [h\ \mathrm{cMpc^{-1}}$]')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('$P_\mathrm{k}$')
    plt.colorbar(label='$P_\mathrm{k}\ [\mathrm{K^2}\ h^{-1}\ \mathrm{cMpc}^3$]')
    plt.savefig(outpath+'ps2d_{}_'.format(kind)+freq_inteval+'.pdf',bbox_inches='tight')
    
    plt.figure()
    plt.pcolormesh(k_per,k_par,ps_err,norm='log',cmap='jet')
    plt.xlim(k_per.min(),k_per.max())
    plt.ylim(k_par.min(),k_par.max())
    plt.xlabel('$k_\perp\ [h\ \mathrm{cMpc^{-1}}$]')
    plt.ylabel('$k_\parallel\ [h\ \mathrm{cMpc^{-1}}$]')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('$\Delta P_\mathrm{k}$')
    plt.colorbar(label='$\Delta P_\mathrm{k}\ [\mathrm{K^2}\ h^{-1}\ \mathrm{cMpc}^3$]')
    plt.savefig(outpath+'ps2d_err_{}_'.format(kind)+freq_inteval+'.pdf',bbox_inches='tight')
    
    # kbins = np.logspace(np.log10(ps_gen.kmin), np.log10(ps_gen.all_k.max()), 50)
    # ps1d = ps_gen.get_ps3d(kbins, data_cube)
    # ps1d.plot(title='Spherically averaged power spectra',kerr_as_kbins=ERROR_P1D)
    # plt.savefig(outpath+'ps1d_'+'.pdf')
    # ps1d.save_to_txt(outpath+'ps1d_'+'.txt')
    
processes = []
for i in res_cubes:
    p = Process(target=cal_ps, args=(i,))
    p.start()
    processes.append(p)
    
for p in processes:
    p.join()
