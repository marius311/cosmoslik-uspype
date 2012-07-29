from pywmap import pywmap
from numpy import zeros
from cosmoslik.modules import Likelihood
import os


class wmap(Likelihood):
    
    def init(self,p):
        self.use = p.get(('wmap','use'),['TT','TE','EE','BB'])
        ttmin, ttmax = p.get(('wmap','TT.lrange'),(2,1200))
        temin, temax = p.get(('wmap','TE.lrange'),(2,800))
        datadir = p.get(('wmap','data_dir'),None)
        if not datadir: raise Exception('Please specify the WMAP data directory in the parameter file with [wmap]{data_dir=..}')
        elif not os.path.exists(datadir): raise Exception("The WMAP data directory you specified in the parameter file does not exist. '%s'"%datadir)
        pywmap.wmapinit(ttmin,ttmax,temin,temax,os.path.normpath(datadir)+'/')
    
    def get_required_models(self, p):
        return ['cl_%s'%x for x in self.use]
    
    def lnl(self, p, model):
        cltt, clte, clee, clbb = [zeros(1202) for _ in range(4)]
        
        for cl,x in zip([cltt,clte,clee,clbb],['TT','TE','EE','BB']):
            if x in self.use:
                m = model['cl_%s'%x]
                s = slice(0,min(len(m),len(cl)))
                cl[s] = m[s]
    
        liketerms = pywmap.wmaplnlike(cltt=cltt[2:],clte=clte[2:],clee=clee[2:],clbb=clbb[2:])
        return sum(liketerms)
    
    
