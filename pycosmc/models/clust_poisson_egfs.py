from pycosmc.modules import Model
from numpy import *

class clust_poisson_egfs(Model):
    
    def get(self, p, required):
        
        def clust_poisson_egfs(spectra, lmax, **kwargs):
            if spectra != 'cl_TT': 
                return zeros(lmax)
            else:
                return p['Aps'] * (arange(lmax)/3000.)**2 + p['Acl'] * (arange(lmax)/3000.)**0.8 

        clust_poisson_egfs.__reduce_ex__ = lambda _: (_unpicklable,(),None,None,None)

        return {'egfs':clust_poisson_egfs}
    
    
def _unpicklable(): pass 
