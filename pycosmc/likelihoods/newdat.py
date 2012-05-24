from pycosmc.modules import Likelihood
from numpy import *
from itertools import islice
from scipy.linalg import cho_factor, cho_solve
import re


class newdat(Likelihood):
    
    xs = ['TT','EE','BB','EB','TE','TB']
    
    def init(self, p):
        with open(p['newdat','file']) as f:
            name = remove_comments(f.next())
            nxs = map(int,remove_comments(f.next()).split())
            has_calib_uncert, calib, calib_err = remove_comments(f.next()).split()
            has_beam_uncertain, beam, beam_err = remove_comments(f.next()).split()
            ilike = int(remove_comments(f.next()))
            self.bands = {}
            
            for x,nx in zip(self.xs,nxs):
                if nx!=0:
                    if f.next().strip()!=x: raise Exception('Error reading newdat file. Expected bandpowers in order %s'%self.xs)
                    self.bands[x] = array([fromstring(remove_comments(s),sep=' ') for s in islice(f,nx)])
                    for _ in islice(f,nx): pass #ignore correlation matrix
        
            self.cov = cho_factor(array([fromstring(remove_comments(s),sep=' ') for s in islice(f,sum(nxs))]))                
                
                
    def get_required_models(self, p):
        return ['cl_%s'%x for x in self.bands.keys()]
                
                
    def lnl(self, p, model):
        
        if p.get('diagnostic',False):
            from matplotlib.pyplot import ion, draw, plot, errorbar, cla, yscale, ylim
            ion()
            cla()
            ells = [mean(lrange) for lrange in self.bands['TT'][:,[5,6]]]
            errorbar(ells,self.bands['TT'][:,1],yerr=self.bands['TT'][:,[3,2]].T,fmt='.')
            plot(ells,[mean(model['cl_TT'][lmin:lmax+1]) for (lmin,lmax) in self.bands['TT'][:,[5,6]]])
            yscale('log')
            ylim(10,6e3)
            draw()
        
        dcl = hstack([array([mean(model['cl_%s'%x][lmin:lmax+1]) for (lmin,lmax) in self.bands[x][:,[5,6]]])-self.bands[x][:,1] for x in self.xs if x in self.bands])
        return dot(dcl,cho_solve(self.cov,dcl))/2
        
        
def remove_comments(s):
    return re.sub('#.*','',s).strip()



#  !File Format:
#  !name
#  !n_bands_TT n_EE, n_BB, n_EB, n_TE, n_TB
#  !has_calib_uncertain calib(amplitude) calib_err(power)
#  !has_beam_uncertain beam beam_err
#  !ilike (0: Gaussian, 1: all x-factor, 2: specified have x-factor)
#  !loop over {
#  ! band-types
#  ! band info: num obs + - x_factor l_min l_max use_x
#  ! correlation-matrix (ignored)
#  ! }  
#  ! covariance matrix
            