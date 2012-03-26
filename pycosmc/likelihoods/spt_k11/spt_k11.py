from numpy import array, fromstring, loadtxt, dot
from scipy.linalg import cho_factor, cho_solve
import os

eff_fr={'dusty':(150,150),
        'radio':(150,150),
        'tsz':(150,150)}

fluxcut = 50

def init(p):
    global spec, sigma, windows, windowrange, datadir
    
    datadir = os.path.join(os.path.dirname(__file__),'bandpowers')
    
    #Load spectrum as covariance
    with open(os.path.join(datadir,'Spectrum_spt20082009.newdat')) as f:
        while 'TT' not in f.readline(): pass
        spec=array([fromstring(f.readline(),sep=' ')[1] for _ in range(47)])
        sigma=cho_factor(array([fromstring(f.readline(),sep=' ') for _ in range(94)])[47:])
        
    #Load windows
    windows = [loadtxt(os.path.join(datadir,'windows','window_0809','window_%i'%i))[:,1] for i in range(1,48)]
    windowrange = (lambda x: slice(min(x),max(x)+1))(loadtxt(os.path.join(datadir,'windows','window_0809','window_1'))[:,0])

    #This likelihood only needs the TT spectrum
    p['models.calculate'].add('cl_TT')
    
    assert p['lmax']>=windowrange.stop-1, "SPT K11 likelihood needs C_ell's to ell=%i"%(p['spt_k11.windowrange'].stop-1)


def lnl(model,p,derivative=0):
    if derivative!=0: raise NotImplementedError("WMAP derivative not implemented yet.")
    
    #Get CMB + foreground model
    cl = model['cl_TT'] 
    if 'fgs' in model: cl += model['fgs'](eff_fr=eff_fr,fluxcut=fluxcut)
    
    #Apply windows and calculate likelihood
    dcl = spec-array([dot(cl[windowrange],w) for w in windows])
    return dot(dcl,cho_solve(sigma, dcl))/2
    