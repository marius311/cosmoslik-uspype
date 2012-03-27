from pywmap import pywmap
from numpy import zeros

def init(p):
    print 'Initializing WMAP likelihood...'
    pywmap.wmapinit()
    for x in p.get('wmap.use',['TT','TE','EE']): p['models.calculate'].add('cl_%s'%x)

def lnl(model,p,derivative=0):
    if derivative!=0: raise NotImplementedError("WMAP derivative not implemented yet.")
    
    cltt, clte, clee, clbb = [zeros(1202) for _ in range(4)]
    
    use = p.get('wmap.use',['TT','TE','EE'])
    if 'EE' in use or 'TE' in use: assert 'TE' in use and 'EE' in use, "Must use either TE and EE or neither."
    
    for cl,x in zip([cltt,clte,clee,clbb],['TT','TE','EE','BB']):
        if x in use:
            s = slice(*p.get('wmap.%s.lrange'%x,(2,1200)))
            cl[s] = model['cl_%s'%x][s]

    return pywmap.wmaplnlike(cltt=cltt[2:],clte=clte[2:],clee=clee[2:],clbb=clbb[2:])
