from pycosmc.modules import Model
from numpy import arange, loadtxt, hstack, pi, exp, zeros, ndarray
import os

class Egfs(Model):
    """
    
    """
    
    def get_egfs(self, p, *args, **kwargs):
        raise NotImplementedError()
    
    
    def get(self, p, required):
        
        def get_egfs(*args, **kwargs):
            
            comps = self.get_egfs(p, *args, **kwargs)
            
            if 'plot' in kwargs:
                from matplotlib.pyplot import subplot
                ax = kwargs.pop('ax',None) or subplot(111)
                colors = self.get_colors(p)
                for comp in (lambda key: comps if key is True else key)(kwargs.pop('plot')):
                    ax.plot(comps[comp],label=comp,color=colors[comp])
                    
            return sum(comps.values())

        get_egfs.__reduce_ex__ = lambda _: (_dont_pickle,(),None,None,None)
        
        return {'egfs':get_egfs}
    
    
def _dont_pickle(): pass 

