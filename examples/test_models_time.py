

import numpy as np
import time
from fitIR import fitIR
from fitIR import models



z = 6.
log10LIR = 10.
lam = np.arange(1.,5000.,1.)


T = 35.
emissivity = 1.6

true_parameters = {'z':z, 'T':T, 'emissivity':emissivity, 'log10LIR': log10LIR}

ts = time.time()

for i in range(1000): models.greybody(T, emissivity).fnu(lam, z=z)
# for i in range(1000): fitIR.model(lam, true_parameters, mod = 'greybody')

te = time.time()

print te-ts

# ---- show plot of observations


T = 35.
emissivity = 1.6
alpha = 2.0

true_parameters = {'z':z, 'T':T, 'emissivity':emissivity, 'alpha': alpha, 'log10LIR': log10LIR}

ts = time.time()

for i in range(1000): models.Casey12(T, emissivity, alpha).fnu(lam, z=z)
# for i in range(1000): fitIR.model(lam, true_parameters, mod = 'Casey12')

te = time.time()

print te-ts





