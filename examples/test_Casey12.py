

import numpy as np
# import seaborn as sns
# import pandas as pd
import matplotlib.pyplot as plt

import os
import sys

sys.path.insert(0, os.path.abspath('..'))

from fitIR import fitIR

import cPickle as pickle



# ------- generate fake observations

lam = np.array([500., 850., 1000., 1300., 2000.])

lam = np.array([500., 1300.])

# lam = np.arange(100.,2000.,100.)

z = 0.01
T = 30.
emissivity = 1.6
alpha = 2.0
log10LIR = 10.
true_parameters = {'z':z, 'T':T, 'emissivity':emissivity, 'alpha': alpha, 'log10LIR': log10LIR}

# obs = fitIR.fake_observations(lam, true_parameters, SNR = 20, mod = 'Casey12')





# ---- show plot of observations

lam = np.arange(1.,5000.,1.)
plt.plot(np.log10(lam), np.log10(fitIR.model(lam, true_parameters, mod = 'Casey12').value))

print fitIR.model(lam, true_parameters, mod = 'Casey12')

# plt.scatter(np.log10(obs.lam), obs.fluxes)
plt.show()



# 
# 
# # ---- initialise source
# 
# source = fitIR.source(obs, mod = 'Casey12') 
# 
# # ---- update priors
# 
# source.prior_def['z'] = {'type': 'delta', 'limits': z}
# # source.prior_def['emissivity'] = {'type': 'delta', 'limits': 1.6}
# 
# 
# # ---- run fitter
# 
# source.runMCMC()
# 
# # ---- save results
# 
# source.save('source.pck')


