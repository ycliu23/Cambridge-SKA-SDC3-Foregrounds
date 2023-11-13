import numpy as np
from ps_eor import datacube
import GPy
from astropy.io import fits
import astropy.constants as const
import pypolychord
from pypolychord.settings import PolyChordSettings
from pypolychord.priors import UniformPrior

start_freq = 106
end_freq = 121
freqs = np.linspace(start_freq,end_freq,151)
inpath = 'Sub-band/images/{start_freq}_{end_freq}MHz/processed_cube_{start_freq}_{end_freq}MHz_msn.h5'.format(start_freq=start_freq,end_freq=end_freq)
data_cube = datacube.CartDataCube.load(inpath)
data = data_cube.data
recomb_data = np.concatenate([data.real,data.imag],axis=1)

kern_smooth = GPy.kern.RBF(1)
kern_mix = GPy.kern.RBF(1)
kern_21 = GPy.kern.Exponential(1)

kern = kern_smooth + kern_mix + kern_21
model = GPy.models.GPRegression(freqs[:,None], recomb_data, kern)

nDims = len(model.optimizer_array) # number of kernel hyperparameters
nDerived = 0

def likelihood(theta):
    
    model['.*var'][0], model['.*var'][1], model['.*var'][2] = theta[0], theta[2], theta[4]
    model['.*len'][0], model['.*len'][1], model['.*len'][2] = theta[1], theta[3], theta[5]
    model['.*Gaussian_noise'] = theta[6]
    
    LML = model.log_likelihood()

    return LML, 0
    
def prior(cube):
    theta = np.zeros_like(cube)
    
    theta[0] = UniformPrior(1e1,1e3)(cube[0])
    theta[1] = UniformPrior(10,1000)(cube[1])
    theta[2] = UniformPrior(1e-1,1e2)(cube[2])
    theta[3] = UniformPrior(0.1,10)(cube[3])
    theta[4] = UniformPrior(1e-5,1e-1)(cube[4])
    theta[5] = UniformPrior(0.001,1)(cube[5])
    theta[6] = UniformPrior(0.01,10)(cube[6])
    
    return theta
    
settings = PolyChordSettings(nDims, nDerived)
settings.base_dir = 'ns_{lower}_{upper}_MHz/'.format(lower=start_freq,upper=end_freq)
settings.file_root = 'gpr'
settings.read_resume = True
# settings.nlive = 175

model_params = ['\sigma^{2}_{smooth}','l_{smooth}','\sigma^{2}_{mix}','l_{mix}','\sigma^{2}_{21}','l_{21}','\sigma^{2}_{n}']
output = pypolychord.run_polychord(likelihood, nDims, nDerived, settings, prior)
paramnames = [('p%i' % i, model_params[i]) for i in range(nDims)]
output.make_paramnames_files(paramnames)
