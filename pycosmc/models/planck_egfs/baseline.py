from pycosmc.modules import Model
from numpy import arange, ones, loadtxt, vstack, hstack, pi
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
            norm_fr = 150
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
        def todl(cl,norm=3000): return (lambda dl: dl/dl[norm])((lambda l: cl*l*(l+1)/2/pi)(arange(cl.shape[0])))
        self.clustered_template =  todl(hstack([[0,0],loadtxt(os.path.join(self.dir,"clustered_150.dat"))[:,1]]),norm=1500)
        self.tsz_template = todl(hstack([[0,0],loadtxt(os.path.join(self.dir,"tsz.dat"))[:,1]]))
        self.ksz_template = todl(hstack([[0,0],loadtxt(os.path.join(self.dir,"ksz_ov.dat"))[:,1]]))
        
    def get(self, p, required):
        
        p_egfs = p.get('egfs',{})
        
        def get_egfs(spectra, fluxcut, freqs, lmax):
            if spectra != 'cl_TT': return 0
            
            fr1, fr2 = freqs
            
            clustered_template = (lambda tilt: hstack([self.clustered_template[:1500]*(1500/3000.)**tilt,(arange(1500,10001)/3000.)**tilt]))(p_egfs['dgcl','tilt'])
            
            return sum([p_egfs['dgpo','amp'] * (arange(lmax)/3000.)**2 * (fr1*fr2/p_egfs['dgpo','norm_fr']**2)**p_egfs['dgpo','alpha'],
                        p_egfs['radio','amp'] * (arange(lmax)/3000.)**2 * (fr1*fr2/p_egfs['radio','norm_fr']**2)**p_egfs['radio','alpha'],
                        p_egfs['dgcl','amp'] * clustered_template[:lmax] * (fr1*fr2/p_egfs['dgcl','norm_fr']**2)**p_egfs['dgcl','alpha'],
                        p_egfs['tsz','amp'] * self.tsz_template[:lmax],
                        p_egfs['ksz','amp'] * self.ksz_template[:lmax]])

        return {'egfs':get_egfs}