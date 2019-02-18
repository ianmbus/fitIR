

import numpy as np
# import seaborn as sns
import matplotlib.pyplot as plt
import cPickle as pickle
import corner
from fitIR import fitIR


source = np.load(open('source.pck','r'))




# ------- calculate medians

samples = source.chain[:, 200:, ].reshape((-1, source.ndim))

for ip, parameter in enumerate(source.parameters):

    print parameter, np.median(samples[:,ip])



# ------- walker figure

walker_fig = True
if walker_fig:

    fig, axarr = plt.subplots(len(source.parameters), sharex=True,figsize=(15,10))
    
    for i in range(0, source.nwalkers):   
        for ip, parameter in enumerate(source.parameters):  
            axarr[ip].plot(source.chain[i,:,ip],alpha=0.3,c='b')

    fig.savefig("walkers.png")



# ------- corner figure

corner_fig = True
if corner_fig:

    # --- determine parameters that are not fixed not fixed

    show = []
    for i, parameter in enumerate(source.parameters):
        if source.prior_def[parameter]['type'] != 'delta': show.append(i)


    if source.obs.truth: 
    
        truths = [source.obs.truth[source.parameters[i]] for i in show]
        fig = corner.corner(samples[:, show], labels = [source.parameters[i] for i in show], truths = truths) #
    
    else:
        
        fig = corner.corner(samples[:, show], labels = [source.parameters[i] for i in show]) #
    
    fig.savefig("triangle.png")



# ------- model figure

model_figure = True
if model_figure:

    color = 'cornflowerblue'

    fig, ax = plt.subplots(nrows=1, sharex=True)

    for params  in samples[np.random.randint(len(samples), size=100)]:
                
        p = {parameter: params[i] for i, parameter in enumerate(source.parameters)} 
               
        lam = np.arange(5.,2000.,1.)
        
        m = fitIR.model(lam, p, source.mod)
        
        ax.plot(lam, m, color="k", alpha=0.02)
    
    
    if source.obs.truth: ax.plot(lam, fitIR.model(lam, source.obs.truth, source.mod), color = color, linewidth=3)
    
    ax.errorbar(source.obs.lam, source.obs.fluxes, yerr=source.obs.errors, fmt='o')
    # ax.scatter(source.obs.lam, source.obs.fluxes, color = color)
    
    
    fig.savefig('model_figure.png')