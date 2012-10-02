from numpy import zeros
from cosmoslik.plugins import Likelihood, SubprocessExtension
import os

class planck_lowell_br(Likelihood):
    """
    
    The Planck low-ell likelihood from Hans-Kristian. Based off of 
    the WMAP likelihood code.

 
    Parameters
    ==========
    
    This module reads the following parameters from the ini file:
    
    [planck_lowell_br].data_dir
    ---------------
        The path to the wmap/data directory.
    
    [planck_lowell_br].gibbs_file
        The path to the Gibbs sigma-l file.

    [planck_lowell_br].use
    ----------
        A subset of ``['TT','TE','EE','BB']`` corresponding to 
        which likelihood terms to use.
        
    [planck_lowell_br].TT.lrange
    ----------------
        The TT range in ell to use in the likelihood
    
    [planck_lowell_br].TE.lrange
    ----------------
        The TE range in ell to use in the likelihood
    
    """
    

    def init(self,p):
        self.use = p.get(('planck_lowell_br','use'),['TT','TE','EE','BB'])
        ttmin, ttmax = p.get(('planck_lowell_br','TT.lrange'),(2,1200))
        temin, temax = p.get(('planck_lowell_br','TE.lrange'),(2,800))
        datadir = p.get(('planck_lowell_br','data_dir'),None)
        gibbsfile = p.get(('planck_lowell_br','gibbs_file'),None)
        if not datadir: raise Exception('Please specify the WMAP data directory in the parameter file with:\n[planck_lowell_br]{\n  data_dir=/path/to/wmap/data\n}')
        elif not os.path.exists(datadir): raise Exception("The WMAP data directory you specified does not exist: '%s'"%datadir)


        if not gibbsfile: raise Exception('Please specify the Gibbs file in the parameter file with:\n[planck_lowell_br]{\n  gibbs_file=/path/to/file\n}')
        elif not os.path.exists(gibbsfile): raise Exception("The Gibbs file you specified does not exist: '%s'"%datadir)

        self.pywmap = SubprocessExtension('pywmap',globals())
        self.pywmap.wmapinit(ttmin,ttmax,temin,temax,os.path.normpath(datadir)+'/',os.path.normpath(gibbsfile))
    
    def get_required_models(self, p):
        return ['cl_%s'%x for x in self.use]
    
    def lnl(self, p, model):
        cltt, clte, clee, clbb = [zeros(1202) for _ in range(4)]
        
        for cl,x in zip([cltt,clte,clee,clbb],['TT','TE','EE','BB']):
            if x in self.use:
                m = model['cl_%s'%x]
                s = slice(0,min(len(m),len(cl)))
                cl[s] = m[s]
    
        liketerms = self.pywmap.wmaplnlike(cltt=cltt[2:],clte=clte[2:],clee=clee[2:],clbb=clbb[2:])
        return sum(liketerms)
    
