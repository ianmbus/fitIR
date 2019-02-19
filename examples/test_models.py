

import numpy as np
import matplotlib.pyplot as plt



import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fitIR.models import *


# This example just plots a pair of models to show the difference between the default Greybody and the Casey12 parameterisation. 
# In both cases the total IR luminosity, temperature, and emissivity should be identical.



lam = np.arange(1.,5000.,1.) # um

# ---- default model = greybody


plt.plot(np.log10(lam), greybody({'T': 40, 'emissivity': 1.6}).Lnu(lam, normalised = True), label = 'greybody')

# ---- Casey12 model (with same T, emissivity)

plt.plot(np.log10(lam),  Casey12({'T': 40, 'emissivity': 1.6, 'alpha': 2.0}).Lnu(lam, normalised = True), label = 'Casey12')


plt.xlabel(r'$\log_{10}(\lambda/\mu m)$')
plt.ylabel(r'$L_{\nu}/erg\, s^{-1}\, Hz^{-1}$')
plt.legend()

plt.show()






# ---------- flux plot

from astropy.cosmology import WMAP9 as cosmo

# --- uses HFLS3 approximately

z = 6.3
lamz = np.arange(1.,5000.,1.) # um

# ---- default model = greybody

p = {'T': 40., 'emissivity': 1.6, 'log10LIR': np.log10(3.) + 13. }

plt.plot(np.log10(lamz), greybody(p).fnu(lamz, z, cosmo))

plt.xlabel(r'$\log_{10}\lambda_{obs}/\mu m$')
plt.ylabel(r'$f_{\nu}/mJy$')

plt.show()




