"""

=================================================
Create your own extra-galactic foreground modules
=================================================

To create your own extra-galactic foreground model, create a subclass
of `pycosmc.models.egfs.Egfs()` and override the function `get_egfs` to return a 
dictionary of extra-galactic foreground components. 

Also passed by keyword to the `get_egfs` function is information from the dataset, such as 

- `spectra` : e.g. ``'cl_TT'`` or ``'cl_EE'``
- `freq` : a dictionary for different effective frequencies, e.g. 
  ``{'dust': 153, 'radio': 151, 'tsz':150}``
- `fluxcut` : the fluxcut in mJy
- `lmax` : the necessary maximum l

Here's an example model, ::

    from pycosmc.models.egfs import Egfs
    
    class MyEgfs(Egfs):
    
        def get(self, p, spectra, freq, fluxcut, lmax, **kwargs):
            return {'single_component': p['amp'] * ones(lmax)}

To use this model, place the above code a file called `MyEgfs.py` in 
the folder `pycosmc/models` and then add `MyEgfs` to the `models` key in the parameter file.
 
"""

from pycosmc.modules import Model
from numpy import arange, loadtxt, hstack, pi, exp, zeros, ndarray
import os

class Egfs(Model):
    
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

