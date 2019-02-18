


import numpy as np
import scipy.stats
import sys
import os
import emcee
import models
import cPickle as pickle


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


    
def model(lam, parameters, mod = 'greybody'):

    if mod == 'greybody':

        z = parameters['z']
        T = parameters['T']
        emissivity = parameters['emissivity']
        log10LIR = parameters['log10LIR']
    
        return 10**log10LIR * models.greybody(T = T, emissivity = emissivity).fnu(lam, z=z).to('uJy')  # return fluxes in uJy   

    if mod == 'Greve12':

        z = parameters['z']
        T = parameters['T']
        emissivity = parameters['emissivity']
        log10LIR = parameters['log10LIR']
    
        return 10**log10LIR * models.Greve12(T = T, emissivity = emissivity).fnu(lam, z=z).to('uJy')  # return fluxes in uJy   


    if mod == 'Casey12':

        z = parameters['z']
        T = parameters['T']
        emissivity = parameters['emissivity']
        alpha = parameters['alpha']
        log10LIR = parameters['log10LIR']
    
        return 10**log10LIR * models.Casey12(T = T, emissivity = emissivity, alpha = alpha).fnu(lam, z=z).to('uJy')  # return fluxes in uJy   






def fake_observations(lam, parameters, SNR = 10., mod = 'greybody'):

    fluxes = model(lam, parameters, mod)

    errors = fluxes/SNR

    fluxes += errors * np.random.randn(len(lam))

    obs = observations(lam, fluxes, errors, truth = parameters)

    return obs


class observations():

    # could instead use a dictionary but I prefer this.

    def __init__(self, lam, fluxes, errors, truth = False):

        self.lam = lam
        self.fluxes = fluxes.to('uJy').value
        self.errors = errors.to('uJy').value
        self.truth = truth
    

class output: pass # --- output class




   
    
    
class source():


    def __init__(self, obs, mod = 'greybody'):
     
     
     
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




    def lnlike(self, fluxes_model):
        """log Likelihood function"""
        
        return -0.5*np.sum(((self.obs.fluxes - fluxes_model)/self.obs.errors)**2 - np.log(self.obs.errors**2))


    def lnprob(self, params):
        """Log probability function"""

        p = {parameter:params[i] for i,parameter in enumerate(self.parameters)}
      
        lp = np.sum([self.priors[parameter].logpdf(p[parameter]) for parameter in self.parameters])
  
        if not np.isfinite(lp):
            return -np.inf

        fluxes_model = model(self.obs.lam, p, mod = self.mod).value


        LNPROB =  lp + self.lnlike(fluxes_model)

        if LNPROB != LNPROB: LNPROB = -np.inf

        return LNPROB


    def runMCMC(self, nwalkers = 50, nsamples = 1000):
    
        ndim = len(self.parameters)
    
        self.update_priors()
    
        self.ndim = ndim
        self.nwalkers = nwalkers        
        self.nsamples = nsamples
    
        pos = [ [self.priors[parameter].rvs() for parameter in self.parameters] for i in range(nwalkers)]
        
        self.sampler = emcee.EnsembleSampler(nwalkers, ndim, self.lnprob, args=())
        self.sampler.run_mcmc(pos, nsamples)
        

        
        
        
        
    def save(self, file):
    
    
        out = output()
        
        out.parameters = self.parameters
        out.mod = self.mod
        
        out.obs = self.obs 
        
        
        out.prior_def = self.prior_def
        
        
        samples = self.sampler.chain[:, 500:, ].reshape((-1, self.ndim))

        out.P = {}
        for ip, parameter in enumerate(self.parameters):
            out.P[parameter] = np.array([np.percentile(samples[:,ip], x) for x in range(0,101)])

        out.ndim = self.ndim 
        out.nwalkers = self.nwalkers 
        out.nsamples = self.nsamples
#         out.chain = self.sampler.chain 

        pickle.dump(out, open(file,'wb'))
