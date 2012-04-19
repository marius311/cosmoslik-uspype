from numpy import *
from collections import namedtuple
import pycosmc.mpi as mpi
from emcee.ensemble import EnsembleSampler
from numpy.random import multivariate_normal

def sample(x,lnl,**kwargs):
    nwalkers = kwargs.get('goodman_weare',{}).get('walkers',100)
    nsamp = kwargs.get('samples',10000)
    
    sampler=EnsembleSampler(nwalkers,len(x), lambda x: -lnl(x,**kwargs)[0], pool=namedtuple('pool',['map'])(mpi.mpi_map))
    p0=mpi.mpi_consistent(multivariate_normal(x,kwargs['_cov'],size=nwalkers))

    for pos,lnprob,_ in sampler.sample(p0,iterations=nsamp/nwalkers):
        for cur_x, cur_lnl in zip(pos.copy(),lnprob.copy()):
            if mpi.is_master(): yield cur_x, 1, cur_lnl, kwargs
    

