from numpy import *
from collections import namedtuple
import pycosmc.mpi as mpi
from emcee.ensemble import EnsembleSampler
from numpy.random import multivariate_normal
from itertools import product
import re

def sample(x,lnl,**kwargs):
    if '_cov' not in kwargs: kwargs['_cov'] = initialize_covariance(kwargs)
    nwalkers = kwargs.get('walkers',100)
    nsamp = kwargs.get('samples',10000)
    
    sampler=EnsembleSampler(nwalkers,len(x), lambda x: -lnl(x,**kwargs)[0], pool=namedtuple('pool',['map'])(mpi.mpi_map))
    p0=mpi.mpi_consistent(multivariate_normal(x,kwargs['_cov'],size=nwalkers))

    for pos,lnprob,_ in sampler.sample(p0,iterations=nsamp/nwalkers):
        for cur_x, cur_lnl in zip(pos.copy(),lnprob.copy()):
            if mpi.is_master(): yield cur_x, 1, cur_lnl, kwargs
    
def get_sampled(params): return params['_sampled']

def initialize_covariance(params):
    """Load the sigma, defaulting to diagonal entries from the WIDTH of each parameter."""
    v=params.get("proposal_matrix","")
    if (v==""): prop_names, prop = [], None
    else: 
        with open(v) as f:
            prop_names = re.sub("#","",f.readline()).split()
            prop = genfromtxt(f)
    sigma = diag([params["*"+name][3]**2 for name in get_sampled(params)])
    common = set(get_sampled(params)) & set(prop_names)
    if common: 
        idxs = zip(*(list(product([ps.index(n) for n in common],repeat=2)) for ps in [get_sampled(params),prop_names]))
        for ((i,j),(k,l)) in idxs: sigma[i,j] = prop[k,l]
    return sigma/len(params['_sampled'])*params.get('proposal_scale',2.4)**2