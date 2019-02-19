

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.cosmology import WMAP9 as cosmo


import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import fitIR
import fitIR.models as models
import fitIR.analyse as analyse

import pickle




z = 6.
T = 35.
emissivity = 1.6
log10LIR = 10.

model = 'greybody'

p_true = {'z':z, 'T':T, 'emissivity':emissivity, 'log10LIR': log10LIR, 'model': model}


# ------- generate fake observations

lamz = np.array([500., 850., 1000., 1300., 2000.])

obs = fitIR.fake_observations(lamz, p_true, cosmo, SNR = 20)

analyse.plot_obs(obs)





# ---- initialise source

source = fitIR.source(obs, mod = 'greybody') 

# ---- set priors

# source.prior_def['z'] = {'type': 'delta', 'value': z} # <---- to fix redshift
source.prior_def['z'] = {'type': 'uniform', 'limits': [2, 8]} 
source.prior_def['log10LIR'] = {'type': 'uniform', 'limits': [6.,15.]}   
source.prior_def['T'] = {'type': 'uniform', 'limits': [20.,60.]} 
source.prior_def['emissivity'] = {'type': 'uniform', 'limits': [1.,2.]} 


# ---- run fitter

# output = source.fit() # --- returns a dictionary with the post-burn chain for each parameter
# pickle.dump(output, open('test_fitter.p','wb'))

output = pickle.load(open('test_fitter.p','rb'))

# ---- analyse results


a = analyse.analyser(output)

P = a.P() # returns 5th, 16th, 50th, 84th, 95th percentile 
for p in P.keys(): print(p, P[p])

a.sed_plot()
a.triangle_plot()









