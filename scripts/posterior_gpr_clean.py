import numpy as np
import GPy
from ps_eor import datacube
import h5py
import os
from anesthetic import read_chains
from multiprocessing import Process

freq_boundaries = np.linspace(106, 196, 7, dtype=int)
freq_pair = [(freq_boundaries[i], freq_boundaries[i+1]) for i in range(6)]
os.system('mkdir -p gpr_processed_cube')

inpath = ''
outpath = ''

def gpr_clean(pair):
    lower, upper = pair
    freqs = np.linspace(lower,upper,151)
    file = '{lower}_{upper}MHz/processed_cube_{lower}_{upper}MHz_msn.h5'.format(lower=lower,upper=upper)
    data_cube = datacube.CartDataCube.load(inpath+file)
    data = data_cube.data
    recomb_data = np.concatenate([data.real, data.imag],axis=1)
    
    samples = read_chains('ns_{}_{}_MHz_estimate/gpr'.format(lower,upper))
    var_smooth, l_smooth, var_mix, l_mix, var_21, l_21, n, _, _, _ = samples.mean()
    
    kern_smooth = GPy.kern.RBF(1,variance=var_smooth,lengthscale=l_smooth)
    kern_mix = GPy.kern.RBF(1,variance=var_mix,lengthscale=l_mix)
    kern_21 = GPy.kern.Exponential(1,variance=var_21,lengthscale=l_21)
    kern = kern_smooth + kern_mix + kern_21
    
    gpr_fit = GPy.models.GPRegression(freqs[:,None], recomb_data, kern)
    gpr_fit.Gaussian_noise.variance.constrain_fixed(n)
    fit, fit_err = gpr_fit.predict(freqs[:, None], kern=gpr_fit.kern.parts[0]+gpr_fit.kern.parts[1], include_likelihood=False)
    fit = fit[:,:fit.shape[1]//2] + 1j*fit[:,fit.shape[1]//2:]
    res = data - fit
    
    res_cube = data_cube.new_with_data(data=res)
    res_cube.save(outpath+'gpr_res_{lower}_{upper}_MHz.h5'.format(lower=lower,upper=upper))
    fit_cube = data_cube.new_with_data(data=fit)
    fit_cube.save(outpath+'gpr_fit_{lower}_{upper}_MHz.h5'.format(lower=lower,upper=upper))
    
    np.savetxt(outpath+'gpr_model_params_{lower}_{upper}_MHz.txt'.format(lower=lower,upper=upper),np.c_[gpr_fit.parameter_names(),gpr_fit.param_array],fmt='%s')
    
processes = []
for i in freq_pair:
    p = Process(target=gpr_clean, args=(i,))
    p.start()
    processes.append(p)
    
for p in processes:
    p.join()
    
