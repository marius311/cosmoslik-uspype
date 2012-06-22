from pycosmc.modules import Model
from numpy import arange, loadtxt, hstack, pi, exp, zeros
import os

class baseline(Model):
    """
    
    [egfs]{
        [dgpo]{
            amp = 10 
            alpha = 3
            norm_fr = 150
        }
        [dgcl]{
            amp = 0
            alpha = 3
            tilt = 1
            norm_fr = 150
        }
        [radio]{
            amp = 0
            alpha = 3
            gamma = 1.2
            norm_fr = 150
            norm_fluxcut = 50
        }
        [tsz]{
            amp = 0
            norm_fr = 150
        }
        [ksz]{
            amp = 0
        }
    
    """
    
    def init(self, p): 
        
        self.dir = os.path.dirname(os.path.abspath(__file__))
        padding = zeros(10000)
        def todl(cl,norm=3000): return (lambda dl: dl/dl[norm])((lambda l: cl*l*(l+1)/2/pi)(arange(cl.shape[0])))
        self.clustered_template =  todl(hstack([[0,0],loadtxt(os.path.join(self.dir,"clustered_150.dat"))[:,1],padding]),norm=1500)
        self.tsz_template = todl(hstack([[0,0],loadtxt(os.path.join(self.dir,"tsz.dat"))[:,1],padding]))
        self.ksz_template = todl(hstack([[0,0],loadtxt(os.path.join(self.dir,"ksz_ov.dat"))[:,1],padding]))
        
        
    def get(self, p, required):
        
        p_egfs = p.get('egfs',{})
        
        def get_egfs(spectra, fluxcut, freqs, lmax, **kwargs):
            if spectra != 'cl_TT': return zeros(lmax)
            
            fr1, fr2 = freqs
            
            tilted_clust = (lambda tilt: hstack([self.clustered_template[:1500]*(1500/3000.)**tilt,(arange(1500,20001)/3000.)**tilt]))(p_egfs['dgcl','tilt'])
            
            return sum([p_egfs['dgpo','amp'] * (arange(lmax)/3000.)**2 * plaw_dep(fr1['dust'], fr2['dust'], p_egfs['dgpo','norm_fr'], p_egfs['dgpo','alpha']),
                        p_egfs['dgcl','amp'] * tilted_clust[:lmax] * plaw_dep(fr1['dust'], fr2['dust'], p_egfs['dgcl','norm_fr'], p_egfs['dgcl','alpha']),
                        p_egfs['radio','amp'] * (fluxcut / p_egfs['radio','norm_fluxcut']) ** (2+p_egfs['radio','gamma']) * (arange(lmax)/3000.)**2 * plaw_dep(fr1['radio'], fr2['radio'], p_egfs['radio','norm_fr'], p_egfs['radio','alpha']),
                        p_egfs['tsz','amp'] * self.tsz_template[:lmax] * tszdep(fr1['tsz'],fr2['tsz'],p_egfs['tsz','norm_fr']),
                        p_egfs['ksz','amp'] * self.ksz_template[:lmax]])
            
        get_egfs.__reduce_ex__ = lambda _: (_unpicklable,(),None,None,None)
        
        return {'egfs':get_egfs}
    
    
def dBdT(fr1,fr0):
    """ dB/dT at T_CMB """
    dBdT,dBdT0 = map((lambda fr: (lambda x0: x0**4 * exp(x0) / (exp(x0)-1)**2)(fr/57.78)),[fr1,fr0])
    return dBdT/dBdT0  
  
def tszdep(fr1,fr2,fr0):
    """The tSZ frequency dependence."""
    t1,t2,t0 = map(lambda fr: (lambda x0: x0*(exp(x0)+1)/(exp(x0)-1) - 4)(fr/56.78),[fr1,fr2,fr0])
    return t1*t2/t0**2

def plaw_dep(fr1,fr2,fr0,alpha):
    """A power-law frequency dependence."""
    return (fr1*fr2/fr0**2.)**alpha / dBdT(fr1,fr0) / dBdT(fr2,fr0)

def _unpicklable(): pass 

