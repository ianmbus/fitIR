

import numpy as np
# import seaborn as sns
# import pandas as pd
import matplotlib.pyplot as plt

from fitIR import fitIR

import cPickle as pickle


z = 6.
log10LIR = 10.
lam = np.arange(1.,5000.,1.)


T = 35.
emissivity = 1.6

true_parameters = {'z':z, 'T':T, 'emissivity':emissivity, 'log10LIR': log10LIR}

ts = time.time()

plt.plot(np.log10(lam), fitIR.model(lam, true_parameters).value)

te = time.time()

print te-ts

# ---- show plot of observations


T = 35.
emissivity = 1.6
alpha = 2.0

true_parameters = {'z':z, 'T':T, 'emissivity':emissivity, 'alpha': alpha, 'log10LIR': log10LIR}

ts = time.time()

plt.plot(np.log10(lam), fitIR.model(lam, true_parameters, model = 'Casey12').value)

te = time.time()

print te-ts

plt.show()



