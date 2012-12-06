from cosmoslik.plugins import Likelihood
from numpy import hstack, arange, array, cumsum, frombuffer, dot, ix_
from scipy.linalg import cho_factor, cho_solve, inv
import cPickle

class camspec_slik(Likelihood):
    
    def init(self,p):
        
        self.labels = [(100,100),(143,143),(217,217),(143,217)]        
        self.freqs = {x:x for x in self.labels}
        self.in_lrange = [(50,1201),(50,2001),(500,2501),(500,2501)]
        self.nl = [u-l for l,u in self.in_lrange]
        self.ells = hstack([arange(*r) for r in self.in_lrange])
        self.out_lrange = p.get(('camspec','lrange'),self.in_lrange)
        
        for ri,ro in zip(self.in_lrange, self.out_lrange):
            if ro and (ro[0]<ri[0] or ro[1]>ri[1]): 
                raise Exception("Camspec lrange outside of available range.")
        
        self.slice = hstack([arange(*(array(r)+s-l)) 
                             for s,(l,u),r in zip([0]+list(cumsum(self.nl)[:-1]),
                                                  self.in_lrange,
                                                  self.out_lrange)
                             if r])
        
        with open(p['camspec','like_file'],'r') as f: self.x, cv = cPickle.load(f)
        self.x *= (self.ells*(self.ells+1))
        cv = ((cv*self.ells*(self.ells+1)).T*self.ells*(self.ells+1)).T
        self.cho_cov = inv(cv[ix_(self.slice,self.slice)])
        
        
    def get_required_models(self, model):
        return ['cl_TT','egfs']
        
    def lnl(self, p, model):
        dcl = self.x[self.slice] - hstack(self.get_cl_model(p, model))[self.slice]
        return dot(dcl,dot(self.cho_cov,dcl))/2
    
    def get_cl_model(self, p, model):
        return [model['cl_TT'][slice(*r)] + 
                model['egfs']('cl_TT',lmax=r[1],freqs=self.freqs[l])[slice(*r)] 
                for r,l in zip(self.in_lrange, self.labels)]
            
#    def plot(self, p, ax=None):
#        from matplotlib.pyplot import figure
#        if ax is None: ax=figure().add_subplot(111)
#        
#        ax.errorbar(arange(*self.lrange),
#                    self.spec[self.lslice],
#                    yerr=self.errorbars[self.lslice])
#        
#        ax.plot(arange(*self.lrange), self.get_cl_model(p, p['_model']))
#        
        
