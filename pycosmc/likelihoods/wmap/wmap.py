from pywmap import pywmap
from numpy import zeros

def init(p):
    print 'Initializing WMAP likelihood...'
    pywmap.wmapinit()
    for x in p.get('wmap.use',['TT','TE','EE']): p['models.calculate'].add('cl_%s'%x)

def lnl(model,p,derivative=0):
    if derivative!=0: raise NotImplementedError("WMAP derivative not implemented yet.")
    
    cltt, clte, clee, clbb = [zeros(1202) for _ in range(4)]
    
    if 'EE' in p['wmap.use']: assert 'TE' in p['wmap.use'], "Can't use EE without TE."
    
    for cl,x in zip([cltt,clte,clee,clbb],['TT','TE','EE','BB']):
        if x in p['wmap.use']:
            s = slice(*p['wmap.%s.lrange'%x])
            cl[s] = model['cl_%s'%x][s]

    return pywmap.wmaplnlike(cltt[2:],clte[2:],clee[2:],clbb[2:])
