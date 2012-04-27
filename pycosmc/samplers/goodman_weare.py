from numpy import *
from collections import namedtuple
import pycosmc.mpi as mpi
from emcee.ensemble import EnsembleSampler
from numpy.random import multivariate_normal
from pycosmc.modules import Sampler


class goodman_weare(Sampler):

    def sample(self,x,lnl,p):
        nwalkers = p.get(('goodman_weare','walkers'),max(100,2*len(x)))
        if nwalkers < 2*len(x): raise Exception("You need at least %i Goodman-Weare walkers (twice the number of sampled parameters)"%(2*len(x)))
        nsamp = p.get('samples',10000)
        
        sampler=EnsembleSampler(nwalkers,len(x), lambda x: -lnl(x,p)[0], pool=namedtuple('pool',['map'])(mpi.mpi_map))
        p0=mpi.mpi_consistent(multivariate_normal(x,p['_cov'],size=nwalkers))
    
        for pos,lnprob,_ in sampler.sample(p0,iterations=nsamp/nwalkers):
            for cur_x, cur_lnl in zip(pos.copy(),lnprob.copy()):
                if mpi.is_master(): yield cur_x, 1, -cur_lnl, p
    

