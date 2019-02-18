

import numpy as np
# import seaborn as sns
# import pandas as pd
import matplotlib.pyplot as plt

from fitIR import fitIR

import cPickle as pickle


# ------- generate fake observations

lam = np.array([500., 850., 1000., 1300., 2000.])

# lam = np.arange(100.,2000.,100.)

z = 6.
T = 35.
emissivity = 1.6
log10LIR = 10.
true_parameters = {'z':z, 'T':T, 'emissivity':emissivity, 'log10LIR': log10LIR}

obs = fitIR.fake_observations(lam, true_parameters, SNR = 20)





# ---- show plot of observations

lam = np.arange(1.,5000.,1.)
plt.plot(np.log10(lam), fitIR.model(lam, true_parameters))

plt.scatter(np.log10(obs.lam), obs.fluxes)
plt.show()





# ---- initialise source

source = fitIR.source(obs, model = 'greybody') 

# ---- update priors

source.prior_def['z'] = {'type': 'delta', 'limits': z}
# source.prior_def['emissivity'] = {'type': 'delta', 'limits': 1.6}


# ---- run fitter

source.runMCMC()

# ---- save results

source.save('source.pck')


