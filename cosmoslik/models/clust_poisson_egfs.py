from cosmoslik.models.egfs import Egfs
from numpy import *

class clust_poisson_egfs(Egfs):
    
    def get_colors(self, p):
        return {'ps':'g','cl':'b'}
        
    def get_egfs(self, p, spectra, lmax, **kwargs):
        if spectra != 'cl_TT': return {}
        else:
            return {'ps': p['Aps'] * (arange(lmax)/3000.)**2,
                    'cl': p['Acl'] * (arange(lmax)/3000.)**0.8}


