
import numpy as np
import scipy.integrate as integrate


class constants: pass # --- constants class

constants.h = 6.626E-34 # J s
constants.c = 2.997E8 # m s^-2
constants.k = 1.380E-23 # J K^-1
constants.L_sol = 3.826E26 * 1E7 # --- luminosity of the Sun in erg s^-1


# ---- build list of redshifts vs. luminosity distances



# redshifts = np.arange(0.,10.,0.01)
# luminosity_distances = np.array([cosmo.luminosity_distance(z).to('cm').value for z in redshifts])
# 
# 





class planck():

    def __init__(self,T):
        
        self.T = T
    
    def Lnu(self, nu):
    
        return (2.*constants.h*(nu**3)*(constants.c**-2))*(1./(np.exp(constants.h*nu/(constants.k*self.T))-1.))


class base():

    def nu(self, lam):
        # --- convert wavelength in um to frequency in Hz
        return constants.c/(lam*1E-6)  
        
    def Lnu(self, lam, normalised = False):
        
        nu = self.nu(lam)
        
        if normalised: 
        
            return self.lnu(nu)/self.normalisation  # --- erg s^-1 Hz^-1, but should be normalised to 1
        
        else:
        
            return constants.L_sol * 10**(self.log10LIR) * self.lnu(nu)/self.normalisation  # --- erg s^-1 Hz^-1
        
  
    def fnu(self, lamz, z, cosmo):   

        # --- return observed frame quantities. THIS IS **NOT** TO BE USED BY THE FITTER DIRECTLY AS IT WOULD BE VERY INEFFICIENT.
    
        luminosity_distance = cosmo.luminosity_distance(z).to('cm').value
    
        lam = lamz/(1. + z)
    
        return  1E3 * 1E23 * (1.+z)*self.Lnu(lam)/(4.*np.pi*luminosity_distance**2) # mJy

   




class greybody(base):

    def __init__(self, p):
    
        if 'log10LIR' not in p:
            p['log10LIR'] = 0.0
    
        self.log10LIR = p['log10LIR']
        self.T = p['T']
        self.emissivity = p['emissivity']
        
        self.kappa = lambda nu: nu**self.emissivity
                
        self.normalisation = integrate.quad(self.lnu, self.nu(1000.), self.nu(8.), full_output=False, limit = 100)[0]
        # self.normalisation /= constants.L_sol.to('erg s^-1').value # --- normalise to solar luminosity


    def lnu(self, nu): # unnormalised lnu
        
        return self.kappa(nu)*planck(self.T).Lnu(nu) 
    

  





class Casey12(base):

    def __init__(self, p):
     
        
        if 'log10LIR' not in p:
            p['log10LIR'] = 0.0
        
        self.log10LIR = p['log10LIR']
        self.T = p['T']
        self.emissivity = p['emissivity']
        self.alpha = p['alpha']
    
        self.N_bb = 1.0
    
        self.lam_0 = 200. #\mu m

        b1 = 26.68
        b2 = 6.246
        b3 = 0.0001905
        b4 = 0.00007243
        L = ( (b1+b2*self.alpha)**-2 + (b3+b4*self.alpha)*self.T)**-1
        self.lam_c = (3./4.)*L
        self.lam_c_m = self.lam_c*10**-6
        
        A1 = (constants.c/self.lam_c_m)**(self.emissivity+3.)
        A2 = np.exp((constants.h*constants.c)/(self.lam_c_m*constants.k*self.T)) - 1.0
        B = self.lam_c**self.alpha
        
        self.N_pl = A1/(A2*B)
        
        self.normalisation = integrate.quad(self.lnu, self.nu(1000.), self.nu(8.), full_output=False, limit = 100)[0]
       
       
    def lnu_lam(self, l):

        l_m = l*10**-6

        PL = self.N_pl*(l**self.alpha)*np.exp(-(l/self.lam_c)**2)

        numerator = (constants.c/l_m)**(self.emissivity+3) 
        denominator = np.exp((constants.h*constants.c)/(l_m*constants.k*self.T)) - 1.0

        BB = self.N_bb*numerator/denominator
        
        return BB + PL  
      
    def lnu(self, nu):
    
        lam = 1E6*constants.c/nu # um
        
        return self.lnu_lam(lam)
    
        

    
        

        
              
















        

# 
# 
# class Greve12():
# 
# 
#     def __init__(self, T, emissivity, nu_c = 3E12):
#     
#         self.T = T
#         self.emissivity = emissivity
#         self.nu_c = nu_c
#                 
#         self.normalisation = integrate.quad(self.Lnunu, self.nu(1000.), self.nu(4.), full_output=False, limit = 100)[0]
#         self.normalisation /= constants.L_sol.to('erg s^-1').value # --- normalise to solar luminosity
#                 
#     def nu(self, lam):
#         # --- convert wavelength in um to frequency in Hz
#         return constants.c/(lam*1E-6)
#  
#  
#           
#     def Lnunu(self, nu):
#         
#         return (1.-np.exp(-(nu/self.nu_c)**self.emissivity))*planck(self.T).Lnu(nu) 
# 
#       
#     def Lnu(self, lam):
#       
#         lam_m = lam*1E-6
#         
#         nu = constants.c/lam_m
# 
#         return self.Lnunu(nu)/self.normalisation * u.erg / (u.s * u.Hz)
# 
#            
#     def fnu(self, lamz, z):
#         
#         # luminosity_distance = cosmo.luminosity_distance(z).to('cm')
#         luminosity_distance = np.interp(z, redshifts, luminosity_distances) * u.cm
#         
#         lam = lamz / (1.+z) # ---- rest-frame wavelength in um
#         
#         fnu = (1.+z)*self.Lnu(lam)/(4.*math.pi*luminosity_distance**2)
#         
#         return fnu
# 
# 


