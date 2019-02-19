

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





# ------- HFLS3 observations from Riechers+

lamz = np.array([250., 350., 500., 880., 1110.])
fluxes = np.array([12., 32.4, 47.3, 33., 21.3])
errors = np.array([2.3, 2.3, 2.8, 2.4, 1.1])
z = 6.34

obs = fitIR.observations(lamz, fluxes, errors, cosmo)

analyse.plot_obs(obs)

# ---- initialise source

source = fitIR.source(obs, mod = 'greybody') 

# ---- set priors

source.prior_def['z'] = {'type': 'delta', 'value': z} # <---- to fix redshift
# source.prior_def['z'] = {'type': 'uniform', 'limits': [2, 8]} 
source.prior_def['log10LIR'] = {'type': 'uniform', 'limits': [6.,15.]}   
source.prior_def['T'] = {'type': 'uniform', 'limits': [20.,60.]} 
source.prior_def['emissivity'] = {'type': 'uniform', 'limits': [1.,2.]} 


# ---- run fitter and save output

# output = source.fit() # --- returns a dictionary with the post-burn chain for each parameter
# pickle.dump(output, open('test_HFLS3.p','wb'))

# ---- analyse results

output = pickle.load(open('test_HFLS3.p','rb')) # --- read in output
# output = pickle.load(open('test_HFLS3_freez.p','rb')) # --- read in output

a = analyse.analyser(output)

P = a.P() # returns 5th, 16th, 50th, 84th, 95th percentile 
for p in P.keys(): print(p, P[p])

a.sed_plot()
a.triangle_plot()









