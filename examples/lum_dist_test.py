

import time
from astropy.cosmology import WMAP9 as cosmo
import numpy as np

ts = time.time()

for i in range(1000): cosmo.luminosity_distance(2.0).to('cm')

te = time.time()

print te - ts


redshifts = np.arange(0.,10.,0.01)
luminosit_distance = np.array([cosmo.luminosity_distance(z).to('cm').value for z in redshifts])
    



ts = time.time()

for i in range(1000): np.interp(2.0432, redshifts, luminosit_distance)

te = time.time()

print te - ts