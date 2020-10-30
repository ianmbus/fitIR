
import numpy as np
import scipy.stats
import sys
import os
import emcee
import matplotlib.pyplot as plt
from . import models

class delta():
    
    def __init__(self, value):     
        self.value = value
        
    def rvs(self):    
        return self.value
        
    def logpdf(self, v):  
    
        if v == self.value:
            return 0.0
        else:
            return -np.inf
        
class custom_pdf():
    
    def __init__(self, cdf_x, cdf_y):
        self.cdf_x = cdf_x
        self.cdf_y = cdf_y
        
    def rvs(self):
        v = np.random.rand(1)
        #can consider rounding this number to a ceratin SF level if code runs to slow
        return (np.interp(v,self.cdf_x,self.cdf_y)[0])
        
    def logpdf(self,v):
        if (v>np.max(self.cdf_y)) or (v<np.min(self.cdf_y)):
            return -np.inf
        else:
            return np.log10(np.interp(v,self.cdf_x,self.cdf_y))


    

def model(lamz, p, cosmo, mod = 'greybody'):

    return getattr(models, mod)(p).fnu(lamz, p['z'], cosmo)  



def fake_observations(lamz, p, cosmo, SNR = 10., mod = 'greybody'):

    fluxes = model(lamz, p, cosmo, mod = mod)

    errors = fluxes/SNR

    fluxes += errors * np.random.randn(len(lamz))

    obs = observations(lamz, fluxes, errors, cosmo, truth = p)

    return obs


class observations():

    def __init__(self, lamz, fluxes, errors, cosmo, truth = False):

        self.cosmo = cosmo
        self.lamz = lamz
        self.fluxes = fluxes
        self.errors = errors
        self.truth = truth
    

class output: pass # --- output class



    
class source():


    def __init__(self, obs, mod = 'greybody'):
     
        
        self.cosmo = obs.cosmo
     
        # --- define observations
     
        self.obs = obs
     
        # --- define model parameters
     
        self.mod = mod
     
        if mod == 'greybody': self.parameters = ['T', 'emissivity', 'z', 'log10LIR']
        if mod == 'Greve12': self.parameters = ['T', 'emissivity', 'z', 'log10LIR']
        if mod == 'Casey12': self.parameters = ['T', 'emissivity', 'alpha', 'z', 'log10LIR']
     
        # --- define default priors
           
        self.prior_def = {}
        
        self.prior_def['z'] = {'type': 'uniform', 'limits': [0.,10.]}
        self.prior_def['log10LIR'] = {'type': 'uniform', 'limits': [6.,15.]}   
        
        if mod == 'greybody':
        
            self.prior_def['T'] = {'type': 'uniform', 'limits': [20.,60.]}
            self.prior_def['emissivity'] = {'type': 'uniform', 'limits': [1.,2.]}
            
        if mod == 'Casey12':
        
            self.prior_def['T'] = {'type': 'uniform', 'limits': [20.,60.]}
            self.prior_def['emissivity'] = {'type': 'uniform', 'limits': [1.,2.]}
            self.prior_def['alpha'] = {'type': 'uniform', 'limits': [1.,2.5]} 
        
        if mod == 'Greve12':
        
            self.prior_def['T'] = {'type': 'uniform', 'limits': [20.,60.]}
            self.prior_def['emissivity'] = {'type': 'uniform', 'limits': [1.,2.]}
            
        
        
        
    def update_priors(self):
    
        # --- update priors

        self.priors = {}

        for parameter in self.parameters:
        
            if self.prior_def[parameter]['type'] == 'uniform':
            
                self.priors[parameter] = scipy.stats.uniform(loc = self.prior_def[parameter]['limits'][0], scale = self.prior_def[parameter]['limits'][1] - self.prior_def[parameter]['limits'][0])  
                
                             
            if self.prior_def[parameter]['type'] == 'delta':
    
                self.priors[parameter] = delta(value = self.prior_def[parameter]['value'])
                
                
            if self.prior_def[parameter]['type'] == 'norm': 
            
                self.priors[parameter] = scipy.stats.norm(loc = self.prior_def[parameter]['loc'], scale = self.prior_def[parameter]['scale'])
                
            if self.prior_def[parameter]['type'] == 'custom_pdf': 
            
                self.priors[parameter] = custom_pdf(cdf_x = self.prior_def[parameter]['cdf_x'], cdf_y = self.prior_def[parameter]['cdf_y'])




    def lnlike(self, fluxes_model):
        """log Likelihood function"""
        
        return -0.5*np.sum(((self.obs.fluxes - fluxes_model)/self.obs.errors)**2 - np.log(self.obs.errors**2))


    def lnprob(self, params):
        """Log probability function"""

        p = {parameter:params[i] for i,parameter in enumerate(self.parameters)}
      
        lp = np.sum([self.priors[parameter].logpdf(p[parameter]) for parameter in self.parameters])

        if not np.isfinite(lp):
            return -np.inf

        fluxes_model = model(self.obs.lamz, p, cosmo = self.cosmo, mod = self.mod)

        LNPROB =  lp + self.lnlike(fluxes_model)

        if LNPROB != LNPROB: LNPROB = -np.inf

        return LNPROB


    def fit(self, nwalkers = 50, nsamples = 1000, burn = 200):
    
        ndim = len(self.parameters)
    
        self.update_priors()
        
        #plt.plot(self.priors['z'].cdf_x,self.priors['z'].cdf_y)
        #plt.show()
    
        self.ndim = ndim
        self.nwalkers = nwalkers        
        self.nsamples = nsamples
    
        p0 = [ [self.priors[parameter].rvs() for parameter in self.parameters] for i in range(nwalkers)]

        
        self.sampler = emcee.EnsembleSampler(nwalkers, self.ndim, self.lnprob, args=())
        pos, prob, state = self.sampler.run_mcmc(p0, burn)
        self.sampler.reset()
        self.sampler.run_mcmc(pos, nsamples)
        
    
        chains = self.sampler.chain[:, :, ].reshape((-1, self.ndim))
    
        samples = {p: chains[:,ip] for ip, p in enumerate(self.parameters)}
        
        
        out = output()
        
        out.obs = self.obs
        out.samples = samples
        out.parameters = self.parameters
        out.model = self.mod
        out.prior_def = self.prior_def
        
        return out

        
        
        
        


        
        
        
        
    
   