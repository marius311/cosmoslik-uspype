from numpy import array, fromstring, loadtxt, dot, arange, diag
from scipy.linalg import cho_factor, cho_solve
import os

eff_fr={'dusty':(150,150),
        'radio':(150,150),
        'tsz':(150,150)}

fluxcut = 50

def init(p):
    global spec, sigma, windows, windowrange, datadir, ells
    
    datadir = os.path.join(os.path.dirname(__file__),'bandpowers')
    
    #Load spectrum as covariance
    with open(os.path.join(datadir,'Spectrum_spt20082009.newdat')) as f:
        while 'TT' not in f.readline(): pass
        spec=array([fromstring(f.readline(),sep=' ')[1] for _ in range(47)])
        sigma=cho_factor(array([fromstring(f.readline(),sep=' ') for _ in range(94)])[47:])
        
    #Load windows
    windows = [loadtxt(os.path.join(datadir,'windows','window_0809','window_%i'%i))[:,1] for i in range(1,48)]
    windowrange = (lambda x: slice(min(x),max(x)+1))(loadtxt(os.path.join(datadir,'windows','window_0809','window_1'))[:,0])
    ells = array([dot(arange(10000)[windowrange],w) for w in windows])
    
    #This likelihood only needs the TT spectrum
    p['_models.get'].add('cl_TT')
    
    assert p['lmax']>=windowrange.stop, "SPT K11 likelihood needs C_ell's to ell=%i"%windowrange.stop


def lnl(model,p,derivative=0):
    if derivative!=0: raise NotImplementedError("WMAP derivative not implemented yet.")
    global cl
    
    #Get CMB + foreground model
    cl = model['cl_TT'] 
    cl += p['Aps']*(arange(len(cl))/3000.)**2 #Hacked PS term until egfs module
    #if 'fgs' in model: cl += model['fgs'](eff_fr=eff_fr,fluxcut=fluxcut)
    cl = array([dot(cl[windowrange],w) for w in windows])
    
    #Apply windows and calculate likelihood
    dcl = spec-cl
    return dot(dcl,cho_solve(sigma, dcl))/2

def diagnostic(axes,p):
    axes['cl_TT'].errorbar(ells,spec,yerr=diag(sigma[0]),fmt='.',label='SPT K11')
    axes['cl_TT'].plot(ells,cl)
    