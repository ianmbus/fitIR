


import numpy as np
import scipy.stats
import sys
import os

import pickle
from . import models
import matplotlib.pyplot as plt
import matplotlib.cm as cm














def plot_obs(obs, save_file = False):

    if obs.truth:

        p = obs.truth

        lamz = np.arange(100.,5000.,1.)

        plt.plot(np.log10(lamz), getattr(models, p['model'])(p).fnu(lamz, p['z'], obs.cosmo), zorder = 0, c='k', alpha = 0.2, lw = 3)

    for i, (l, f, e) in enumerate(zip(obs.lamz, obs.fluxes, obs.errors)):

        c = cm.viridis(i/len(obs.lamz))
        plt.scatter(np.log10(l), f, c=[c], zorder = 1)
        plt.plot([np.log10(l)]*2, [f-e, f+e], c=c, zorder = 1)

    plt.xlabel(r'$\lambda_{obs}/\mu m$')
    plt.ylabel(r'$f_{\nu}/mJy$')

    plt.show()







   
        
        
class analyser():


    def __init__(self, output):
    
        self.obs = output.obs
        self.s = output.samples
        self.params = output.parameters
        self.model = output.model
        self.prior_def = output.prior_def

        self.p_fit = self.median()

            
    def median(self):
    
        return {p: np.median(self.s[p]) for p in self.params}   

    def P(self):
    
        return {p: np.array([np.percentile(self.s[p], x) for x in [5., 16., 50., 84., 95.]]) for p in self.params}   
        
    
    def sed_plot(self, save_file = False):
    
        # --- plot full input model

        lamz = np.arange(100.,5000.,1.)

        if self.obs.truth:

            plt.plot(np.log10(lamz), getattr(models, self.model)(self.obs.truth).fnu(lamz, self.obs.truth['z'], self.obs.cosmo), zorder = 0, c='k', alpha = 0.2, lw = 3, label = 'truth')

        plt.plot(np.log10(lamz), getattr(models, self.model)(self.p_fit).fnu(lamz, self.p_fit['z'], self.obs.cosmo), zorder = 0, c='k', alpha = 0.5, lw = 2, ls=':', label = 'fit')


        for i, (l, f, e) in enumerate(zip(self.obs.lamz, self.obs.fluxes, self.obs.errors)):

            c = cm.viridis(i/len(self.obs.lamz))
            plt.scatter(np.log10(l), f, c=[c], zorder = 1)
            plt.plot([np.log10(l)]*2, [f-e, f+e], c=c, zorder = 1)


        plt.xlabel(r'$\lambda_{obs}/\mu m$')
        plt.ylabel(r'$f_{\nu}/mJy$')
        plt.legend()

        if save_file:
            plt.savefig(save_file)

        plt.show()



    def triangle_plot(self, save_file = False, contours = False, hist2d = True, use_prior_range = True, bins = 50):
    
        
        # only plot parameters that actually vary
        
        self.pparams = []
        
        for p in self.params: 
            if self.prior_def[p]['type']=='uniform': self.pparams.append(p)
    
        n = len(self.pparams)
    
        # ---- initialise figure

        plt.rcParams['mathtext.fontset'] = 'stixsans'
        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.size'] = 7 # perhaps should depend on number of parameters to be plotted

        plt.rcParams['ytick.labelsize'] = 4 # perhaps should depend on number of parameters to be plotted
        plt.rcParams['xtick.labelsize'] = 4 # perhaps should depend on number of parameters to be plotted
    
        plt.rcParams['ytick.direction'] = 'in'    # direction: in, out, or inout
        plt.rcParams['xtick.direction'] = 'in'    # direction: in, out, or inout
    
        plt.rcParams['ytick.minor.visible'] = True
        plt.rcParams['xtick.minor.visible'] = True
    

        fig, axes = plt.subplots(n,n, figsize = (8,8))

        left  = 0.125  # the left side of the subplots of the figure
        right = 0.9    # the right side of the subplots of the figure
        bottom = 0.1   # the bottom of the subplots of the figure
        top = 0.9      # the top of the subplots of the figure
        wspace = 0.02   # the amount of width reserved for blank space between subplots
        hspace = 0.02   # the amount of height reserved for white space between subplots
    
        fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)



        # ---- loop over parameters

        for i in np.arange(n):
            for j in np.arange(n):
                
                
                # axes[i,j].text(0.5,0.5,str(i)+str(j), transform=axes[i,j].transAxes) # label panels
                
                axes[i,j].locator_params(axis = 'x', nbins=3)
                axes[i,j].locator_params(axis = 'y', nbins=3)
      
                pi = self.pparams[i]
                pj = self.pparams[j]
                

            
                if i!=0 and j==0 and j<n-1:
                
                    axes[i,j].set_ylabel(pi)
                    
                if i==(n-1):
                
                    axes[i,j].set_xlabel(pj)


                


                if j == i:
                
                    median = np.percentile(self.s[pi], 50)
                
                    if use_prior_range:
                    
                        range = self.prior_def[pi]['limits']
                    
                    else:
                
                        IQR = np.percentile(self.s[pi], 75) - np.percentile(self.s[pi], 25)
                        
                        range = [median-3*IQR, median+3*IQR]

                    
                
                    N,b,p = axes[i,j].hist(self.s[pi], bins = bins, range = range, color = '0.7', edgecolor = '0.7')
                
                    mxN = np.max(N)
                
                    if self.obs.truth: axes[i,j].axvline(self.obs.truth[pj], lw = 1, c = 'k', alpha = 0.2)
                
                    axes[i,j].scatter(median,mxN*1.3,c='k',s=5)
                
                
                    axes[i,j].plot([np.percentile(self.s[pi], 16.), np.percentile(self.s[pi], 84.)],[mxN*1.3]*2,c='k',lw=1)
                
                    axes[i,j].set_xlim(range)
                    axes[i,j].set_ylim([0.,mxN*2.])
                    
                    axes[i,j].spines['right'].set_visible(False)
                    axes[i,j].spines['top'].set_visible(False)
                    axes[i,j].spines['left'].set_visible(False)
                    axes[i,j].yaxis.set_ticks_position('none')
                    axes[i,j].axes.get_yaxis().set_ticks([])


                elif j < i:
                    
                    if self.obs.truth:
                    
                        if hist2d: 
                            c = '1.0'
                        else:
                            c = 'k'
                        
                        axes[i,j].axhline(self.obs.truth[pi], lw = 1, c = c, alpha = 0.2)
                        axes[i,j].axvline(self.obs.truth[pj], lw = 1, c = c, alpha = 0.2)
                    
    
                    if use_prior_range:

                        rangei = self.prior_def[pi]['limits']
                        rangej = self.prior_def[pj]['limits']                    
                    
                    else:
    
                        IQR = np.percentile(self.s[pi], 75) - np.percentile(self.s[pi], 25)
                        median = np.percentile(self.s[pi], 50)
                        rangei = [median-3*IQR, median+3*IQR]
    
                        IQR = np.percentile(self.s[pj], 75) - np.percentile(self.s[pj], 25)
                        median = np.percentile(self.s[pj], 50)
                        rangej = [median-3*IQR, median+3*IQR]
                    
    
                    H, xe, ye = np.histogram2d(self.s[pj], self.s[pi], bins = bins, range = [rangej, rangei]) 
                        
                    H = H.T
  
                    xlims = [xe[0], xe[-1]]
                    ylims = [ye[0], ye[-1]]
    
                    axes[i,j].set_xlim(xlims)
                    axes[i,j].set_ylim(ylims)
    
                    if hist2d: 
                        X, Y = np.meshgrid(xe, ye)
                        axes[i,j].pcolormesh(X, Y, H, cmap = 'plasma',linewidth=0,rasterized=True) #

                    if j != 0: axes[i,j].set_yticklabels([])


                    # --- add contours
                    
                    if contours: 
                    
                        norm=H.sum() # Find the norm of the sum
                        # Set contour levels
                        # contour1=0.99 
                        contour2=0.95
                        contour3=0.68

                        # Set target levels as percentage of norm
                        # target1 = norm*contour1
                        target2 = norm*contour2
                        target3 = norm*contour3

                        # Take histogram bin membership as proportional to Likelihood
                        # This is true when data comes from a Markovian process
                        def objective(limit, target):
                            w = np.where(H>limit)
                            count = H[w]
                            return count.sum() - target

                        # Find levels by summing histogram to objective
                        # level1= scipy.optimize.bisect(objective, H.min(), H.max(), args=(target1,))
                        level2= scipy.optimize.bisect(objective, H.min(), H.max(), args=(target2,))
                        level3= scipy.optimize.bisect(objective, H.min(), H.max(), args=(target3,))

                        # For nice contour shading with seaborn, define top level
                        level4=H.max()
                        levels=[level2,level3]

                        bwi = (xe[1]-xe[0])/2.
                        bwj = (ye[1]-ye[0])/2.

                        if hist2d: ccolor = '1.0'

                        axes[i,j].contour(xe[:-1]+bwi, ye[:-1]+bwj, H, levels=levels, linewidths=0.5, colors=ccolor)
                        
                        
                        

                else:
               
                    axes[i,j].set_axis_off() 



        if save_file: fig.savefig(save_file, dpi = 300)

        plt.show()
        fig.clf()        
        
        
        
        
    
   