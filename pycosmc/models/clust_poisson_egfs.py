from pycosmc.modules import Model
from numpy import *

class clust_poisson_egfs(Model):
    
    def get(self, p, required):
        
        def clust_poisson_egfs(spectra, lmax, **kwargs):
            if spectra != 'cl_TT': 
                return zeros(lmax)
            else:
                return p['Aps'] * (arange(lmax)/3000.)**2 + p['Acl'] * (arange(lmax)/3000.)**0.8 

        return {'egfs':clust_poisson_egfs}