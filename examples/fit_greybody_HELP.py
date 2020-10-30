print('importing modules')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.cosmology import WMAP9 as cosmo
from astropy.table import Table,Column,join
import time

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname('.'), '..')))

import fitIR
import fitIR.models as models
import fitIR.analyse as analyse

import pickle

import h5py
import glob
print('done importing')

t1 = time.time()
taskid = np.int(os.environ['SGE_TASK_ID'])-1
batch_size = 50
low = taskid*batch_size
up = taskid*batch_size + (batch_size)

files_done = glob.glob('../data/*.pkl')
nums_done = []
for file in files_done:
    nums_done.append(int(file.replace('.pkl','').split('_')[-1]))
if taskid in nums_done:
    print('done already')
    sys.exit()
                   

inputs = h5py.File('../../greybody_input_EN1.h5','r')
#inputs = Table.read('../data/greybody_input_EN1.fits')
print('read in table')

if low>len(inputs['help_id']):
    sys.exit()

if up>len(inputs['help_id']):
    up = len(inputs['help_id'])-1
             
             
ids = inputs['help_id'][low:up]
nums = inputs['num'][low:up]

names = ['Temperature_l','Temperature','Temperature_u','redshift_l','redshift','redshift_u','log10lir_l','log10lir','log10lir_u']

f_100 = inputs['f_pacs_green'][low:up]/1E6
f_160 = inputs['f_pacs_red'][low:up]/1E6
f_250 = inputs['f_spire_250'][low:up]/1E3
f_350 = inputs['f_spire_350'][low:up]/1E3
f_500 = inputs['f_spire_500'][low:up]/1E3

ferr_100 = inputs['ferr_pacs_green'][low:up]/1E6
ferr_160 = inputs['ferr_pacs_red'][low:up]/1E6
ferr_250 = inputs['ferr_spire_250'][low:up]/1E3
ferr_350 = inputs['ferr_spire_350'][low:up]/1E3
ferr_500 = inputs['ferr_spire_500'][low:up]/1E3

lamz = np.array([100., 160., 250., 350., 500.]) # um
         
outputs = []
print('starting fitting')
for n,num in enumerate(nums):
    print(n)
    fluxes = np.array([f_100[n], f_160[n], f_250[n], f_350[n], f_500[n]])
    print(fluxes)
    errors = np.array([ferr_100[n], ferr_160[n], ferr_250[n], ferr_350[n], ferr_500[n]])
    print(errors)
    
    mask = ~np.isnan(fluxes)
    
    obs = fitIR.observations(lamz[mask], fluxes[mask], errors[mask], cosmo)

    source = fitIR.source(obs, mod = 'greybody') 

    #read in the z pdf here
    z_pdf = h5py.File('../../HELP/dmu_products/dmu24/dmu24_ELAIS-N1/data/pz_hb_en1.hdf', 'r')
    cdf_y = np.cumsum(z_pdf['pz'][num])/np.sum(z_pdf['pz'][num])
    cdf_x = z_pdf['zgrid'].value
    z_pdf.close()
    
    if np.sum(np.isnan(cdf_y))==len(cdf_y):
        outputs.append(np.nan)
        print('no redshift pdf')
        continue

    source.prior_def['z'] = {'type': 'custom_pdf', 'cdf_x': cdf_y, 'cdf_y':cdf_x} # <---- to fix redshift
    #source.prior_def['z'] = {'type': 'uniform', 'limits': [2, 8]} 
    source.prior_def['log10LIR'] = {'type': 'uniform', 'limits': [8.,14.]}   
    source.prior_def['T'] = {'type': 'uniform', 'limits': [20.,60.]} 
    source.prior_def['emissivity'] = {'type': 'delta', 'value': 1.5} 

    output = source.fit()
    outputs.append([output,ids[n]])
    '''a = analyse.analyser(output)

    P = a.P()
    
    tmp = []
    for key in P.keys():
         if key=='emissivity':
             continue
         tmp.append(P[key][1])
         tmp.append(P[key][2])
         tmp.append(P[key][3])
    output.append(tmp)

output = np.array(output)
tbl_out = Table()
for n in range(len(output[0])):
    col = Column(name=names[n],data=output[:,n])
    tbl_out.add_column(col)'''
    
file = open('../data/lir_help_{}.pkl'.format(taskid),'wb')
pickle.dump(outputs,file)
file.close()
#NEED TO WRITE THE TABLE ONCE THE FILE STRUCTURE HAS BEEN DECIDED
t2 = time.time()

print('time taken is: {}'.format(t2-t1))
