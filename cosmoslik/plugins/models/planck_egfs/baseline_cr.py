from cosmoslik.plugins import Model
from .. egfs import egfs
from numpy import arange, loadtxt, hstack, pi, exp, zeros, sqrt
import os

class baseline_cr(egfs):
    """
    
    Notes:
    
    *amplitudes are D_ell at 3000
    
    [egfs]{
        [low_fr]{
            tied_dusty_alpha=True
            [dgpo]{
                amp = 6
                alpha = 3.5 
                norm_fr = 143
            }
            [dgcl]{
                amp_lin = 1
                amp_nonlin = 1 
                norm_fr = 143  
            }
            [radio]{
                amp = 60 
                alpha = -0.5 
                gamma = -1 
                norm_fr = 143
                norm_fluxcut = 200
            }
            [tsz]{
                amp = 5
                norm_fr = 143
            }
            [ksz]{
                amp = 0 
            }
        }
    
        [high_fr_353]{
            cor = 0.8
            [dgpo]{
                amp = 1e6
            }
            [dgcl]{
                amp_lin = 1e6 
                amp_nonlin = 1e6
            }
        }
        
        [high_fr_545]{
            cor = 0.6
            [dgpo]{
                amp = 1e8
            }
            [dgcl]{
                amp_lin = 1e8
                amp_nonlin = 1e8
            }
        }

        [high_fr_857]{
            cor = 0.4
            [dgpo]{
                amp = 1e10
            }
            [dgcl]{
                amp_lin = 1e10
                amp_nonlin = 1e10
            }
        }
    }
    
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
        
        lowp = p.get(('egfs','low_fr'),{})
        
        if lowp.get('tied_dusty_alpha',False): lowp['dgcl','alpha'] = lowp['dgpo','alpha']
        
        fr1, fr2 = freqs
        
        def fr2label(fr):
            if 300<fr<500: return '353'
            elif 500<fr<800: return '545'
            else: return '857'
        
        dustcomp = {}
        
        for i,fr in enumerate(freqs):
            if fr['dust']<300:
                dustcomp[i] = {'dgpo': sqrt(lowp['dgpo','amp']) * (arange(lmax)/3000.) * plaw_dep(fr['dust'], lowp['dgpo','norm_fr'], lowp['dgpo','alpha']),
                               'dgcl': sqrt(lowp['dgcl','amp_lin'] * self.clustered_template[:lmax]) * plaw_dep(fr['dust'], lowp['dgcl','norm_fr'], lowp['dgcl','alpha']) + 
                                            lowp['dgcl','amp_nonlin'] * (arange(lmax)/self.norm_ell)**p.get(('dgcl','tilt' ),0.8) * plaw_dep(fr['dust'], lowp['dgcl','norm_fr'], lowp['dgcl','alpha'])}
            else:
                highp = p['egfs','high_fr_%s'%fr2label(fr['dust'])]

                dustcomp[i] = {'dgpo': sqrt(highp['dgpo','amp']) * (arange(lmax)/3000.), 
                               'dgcl': sqrt(highp['dgcl','amp_lin'] * self.clustered_template[:lmax] +
                                            highp['dgcl','amp_nonlin'] * (arange(lmax)/self.norm_ell)**p.get(('dgcl','tilt' ),0.8))}


        ffr1, ffr2 = fr1['dust'], fr2['dust']
        
        if min([ffr1,ffr2])<300 and max([ffr1,ffr2])>300: 
            corr = p['egfs','high_fr_%s'%fr2label(max([ffr1,ffr2])),'cor']
        else: 
            corr=1
        
        comps = {}
        for x in ['dgpo','dgcl']: comps[x] = corr*dustcomp[0][x]*dustcomp[1][x]
        comps.update({'radio': lowp['radio','amp'] * (fluxcut / lowp['radio','norm_fluxcut']) ** (2+lowp['radio','gamma']) * (arange(lmax)/3000.) * plaw_dep2(fr1['radio'], fr2['radio'], lowp['radio','norm_fr'], lowp['radio','alpha']),
                      'tsz': lowp['tsz','amp'] * self.tsz_template[:lmax] * tszdep(fr1['tsz'],fr2['tsz'],lowp['tsz','norm_fr']),
                      'ksz': lowp['ksz','amp'] * self.ksz_template[:lmax]})
            
        return comps
    
    
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

