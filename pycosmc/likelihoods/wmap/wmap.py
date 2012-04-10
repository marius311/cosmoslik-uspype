from pywmap import pywmap
from numpy import zeros, loadtxt
import os


def init(p):
    global ttdat, use
    print 'Initializing WMAP likelihood...'
    
    use = p.get('wmap.use',['TT','TE','EE','BB'])
    if 'EE' in use or 'TE' in use: assert 'TE' in use and 'EE' in use, "Must use either TE and EE, or neither."
    for x in use: p['_models.get'].add('cl_%s'%x)
    ttmin, ttmax = p.get('wmap.TT.lrange',(2,1200))
    temin, temax = p.get('wmap.TE.lrange',(2,800))
    pywmap.wmapinit(ttmin,ttmax,temin,temax)


def lnl(model,p,derivative=0):
    if derivative!=0: raise NotImplementedError("WMAP derivative not implemented yet.")
    
    cltt, clte, clee, clbb = [zeros(1202) for _ in range(4)]
    
    for cl,x in zip([cltt,clte,clee,clbb],['TT','TE','EE','BB']):
        if x in use:
            m = model['cl_%s'%x]
            s = slice(0,min(len(m),len(cl)))
            cl[s] = m[s]

    liketerms = pywmap.wmaplnlike(cltt=cltt[2:],clte=clte[2:],clee=clee[2:],clbb=clbb[2:])
    return sum(liketerms)


#def diagnostic(axes,p):
##    axes['cl_TT'].errorbar(ttdat[:,0],ttdat[:,3],yerr=ttdat[:,4],fmt='.',label='WMAP7')
#    axes['cl_TT'].plot(p['_model']['cl_TT'][:1200])
