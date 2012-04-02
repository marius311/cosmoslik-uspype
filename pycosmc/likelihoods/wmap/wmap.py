from pywmap import pywmap
from numpy import zeros, loadtxt
import os

def init(p):
    global ttdat, use
    print 'Initializing WMAP likelihood...'
    
    use = p.get('wmap.use',['TT','TE','EE'])
    if 'EE' in use or 'TE' in use: assert 'TE' in use and 'EE' in use, "Must use either TE and EE, or neither."
    
    pywmap.wmapinit()
    for x in p.get('wmap.use',['TT','TE','EE']): p['_models.get'].add('cl_%s'%x)
    
    ttdat = loadtxt(os.path.join("/home/marius/workspace/mspec/dat/external/wmap_binned_tt_spectrum_7yr_v4p1.txt"))


def lnl(model,p,derivative=0):
    if derivative!=0: raise NotImplementedError("WMAP derivative not implemented yet.")
    
    cltt, clte, clee, clbb = [zeros(1202) for _ in range(4)]
    
    for cl,x in zip([cltt,clte,clee,clbb],['TT','TE','EE','BB']):
        if x in use:
            s = slice(*p.get('wmap.%s.lrange'%x,(2,1200)))
            cl[s] = model['cl_%s'%x][s]

    return pywmap.wmaplnlike(cltt=cltt[2:],clte=clte[2:],clee=clee[2:],clbb=clbb[2:])


def diagnostic(axes,p):
    axes['cl_TT'].errorbar(ttdat[:,0],ttdat[:,3],yerr=ttdat[:,4],fmt='.',label='WMAP7')
    axes['cl_TT'].plot(p['_model']['cl_TT'][:1200])
