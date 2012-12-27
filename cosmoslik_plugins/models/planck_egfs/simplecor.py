from cosmoslik.plugins.models.egfs import egfs
from numpy import arange, hstack, loadtxt, zeros, exp, sqrt, pi, mean
import os


from cosmoslik.plugins import Model
from .. egfs import egfs
from numpy import arange, loadtxt, hstack, pi, exp, zeros, sqrt
import os

class simplecor(egfs):
    """
    
    """
    
    def init(self, p): 
        self.dir = os.path.dirname(os.path.abspath(__file__))
        padding = zeros(10000)
        self.norm_ell = float(p.get('egfs',{}).get('norm_ell',3000))
        def todl(cl,norm=self.norm_ell): return (lambda dl: dl/dl[norm])((lambda l: cl*l*(l+1)/2/pi)(arange(cl.shape[0])))
        self.clustered_template =  todl(hstack([[0,0],loadtxt(os.path.join(self.dir,"clustered_150.dat"))[:,1],padding]))
        self.tsz_template = todl(hstack([[0,0],loadtxt(os.path.join(self.dir,"tsz.dat"))[:,1],padding]))
        self.ksz_template = todl(hstack([[0,0],loadtxt(os.path.join(self.dir,"ksz_ov.dat"))[:,1],padding]))
        
        if p.get(('egfs','tied_dusty_alpha'),False):
            if ('egfs','dgcl.alpha') in p and ('egfs','dgpo.alpha') in p:
                raise Exception("When setting tied_dusty_alpha=True delete egfs.dgcl.alpha") 



    def get_colors(self, p):
        return {'dgpo':'g','dgcl':'g','radio':'orange','tsz':'magenta','ksz':'cyan'}

    def get_egfs(self, p, spectra, fluxcut, freqs, lmax, **kwargs):
        if spectra != 'cl_TT': return zeros(lmax)
        
        freqs = tuple(fr if isinstance(fr,dict) else {'tsz':fr} for fr in freqs)
        
        frlbl = tuple([n for n,(l,u) in {100:(80,120),
                                         143:(130,160), 
                                         217:(200,240)}.items() 
                       if l<mean(fr.values())<u][0] for fr in freqs)

        comps={}
        
        p = p['egfs']
        
        comps['tsz'] = p['a_tsz'] * self.tsz_template[:lmax] * tszdep(freqs[0]['tsz'],freqs[1]['tsz'], p['tsz_norm_fr']) 
        comps['ksz'] = p['a_ksz'] * self.ksz_template[:lmax] 
        
        if frlbl==(100,100):
            comps['ps_100'] = p['a_ps_100']*(arange(lmax)/self.norm_ell)**2
        elif frlbl==(143,143):
            comps['ps_143'] = p['a_ps_143']*(arange(lmax)/self.norm_ell)**2
            comps['cib_143'] = p['a_cib_143']*(arange(lmax)/self.norm_ell)**0.7
            comps['tsz_cib'] = - p['xi'] * sqrt(p['a_tsz'] * self.tsz_template[:lmax] * tszdep(143,143,p['tsz_norm_fr']) * p['a_cib_143'] * (arange(lmax)/self.norm_ell)**0.7)
        elif frlbl==(217,217):
            comps['ps_217'] = p['a_ps_217']*(arange(lmax)/self.norm_ell)**2
            comps['cib_217'] = p['a_cib_217']*(arange(lmax)/self.norm_ell)**0.7
        elif tuple(sorted(frlbl))==(143,217):
            comps['ps_143_217'] = p['r_ps']*sqrt(p['a_ps_143']*p['a_ps_217'])*(arange(lmax)/self.norm_ell)**2
            comps['cib_143_217'] = p['r_cib']*sqrt(p['a_cib_143']*p['a_cib_217'])*(arange(lmax)/self.norm_ell)**0.7
            comps['tsz_cib'] = - p['xi'] * sqrt(p['a_tsz'] * self.tsz_template[:lmax] * tszdep(143,143,p['tsz_norm_fr']) * p['a_cib_217'] * (arange(lmax)/self.norm_ell)**0.7)
            
        return comps

A_radio
alpha_radio

A_dps_143 
A_dps_217
A_dps_545

A_lin_cl_143
A_lin_cl_217
A_lin_cl_545

A_nonlin_cl_143
A_nonlin_cl_217
A_nonlin_cl_545

rcib_143_217
rcib_217_545
rcib_143_545

A_tsz
A_ksz

r_tsz_cib


    
    
    
def dBdT(fr1,fr0):
    """ dB/dT at T_CMB """
    dBdT,dBdT0 = map((lambda fr: (lambda x0: x0**4 * exp(x0) / (exp(x0)-1)**2)(fr/56.78)),[fr1,fr0])
    return dBdT/dBdT0  
  
def tszdep(fr1,fr2,fr0):
    """The tSZ frequency dependence."""
    t1,t2,t0 = map(lambda fr: (lambda x0: x0*(exp(x0)+1)/(exp(x0)-1) - 4)(fr/56.78),[fr1,fr2,fr0])
    return t1*t2/t0**2

def plaw_dep(fr,fr0,alpha):
    """A power-law frequency dependence."""
    return (float(fr)/fr0)**alpha / dBdT(fr,fr0)

def plaw_dep2(fr1,fr2,fr0,alpha):
    """A power-law frequency dependence."""
    return (float(fr1)*fr2/fr0**2)**alpha / dBdT(fr1,fr0) / dBdT(fr2,fr0)


def _unpicklable(): pass 
        
