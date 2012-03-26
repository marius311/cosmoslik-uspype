import emcee, re
from numpy import *
from collections import namedtuple
from mpi import *
from itertools import product


class NamedEnsembleSampler(emcee.EnsembleSampler):
    def __init__(self, nwalkers, params, lnprob, extra_params={},**kwargs):
        self.params = params
        dim = len(params)
        lnprob2 = lambda x,*args: -lnprob(dict(extra_params,**dict(zip(params,x))),*args)
        super(NamedEnsembleSampler, self).__init__(nwalkers,dim,lnprob2,**kwargs)
        
    def sample(self, pos0, **kwargs):
        return super(NamedEnsembleSampler, self).sample(vstack([pos0[k] for k in self.params]).T,**kwargs)


def get_varied(params): return params["$VARIED"]
def get_outputted(params): return params["$OUTPUT"]

def initialize_covariance(params):
    """
    Load the covariance, defaulting to diagonal entries from the WIDTH of each parameter
    """
    
    v=params.get("proposal_matrix","")
    if (v==""): prop_names, prop = [], None
    else: 
        with open(v) as file:
            prop_names = re.sub("#","",file.readline()).split()
            prop = genfromtxt(file)

    params["$COV"] = diag([params["*"+name][3]**2 for name in get_varied(params)])
    common = set(get_varied(params)) & set(prop_names)
    if common: 
        idxs = zip(*(list(product([ps.index(n) for n in common],repeat=2)) for ps in [get_varied(params),prop_names]))
        for ((i,j),(k,l)) in idxs: params["$COV"][i,j] = prop[k,l]

            


def mcmc(p, lnl, mpi=False):
    nwalkers = p.get('walkers',100)
    nsamp = p.get('samples',10000)
    initialize_covariance(p)
    
    lnl2 = lambda p: lnl(p) if all([p['*'+k][1]<p[k]<p['*'+k][2] for k in get_varied(p)]) else inf
    
    sampler=NamedEnsembleSampler(nwalkers,get_varied(p),lnl2,extra_params=p,pool=namedtuple('pool',['map'])(mpi_map) if mpi else None)
    
    p0=mpi_consistent(dict(zip(get_varied(p),random.multivariate_normal([p[k] for k in get_varied(p)],p['$COV'],size=nwalkers).T)))

    if 'chain' in p and is_master(): file=open(p['chain'],'w')
    else: file=None
    
    if file: file.write('#'+' '.join(get_varied(p))+'\n')
    for i,(pos,lnprob,state) in enumerate(sampler.sample(p0,iterations=nsamp/nwalkers),1):
        if file:
            savetxt(file,pos)
            file.flush()
            
        if is_master(): print 'steps=%i approval=%.3f best=%.3f'%(i*nwalkers,sampler.acceptance_fraction.mean(),-sampler.lnprobability[:,:i].max())
    
    if file: file.close()
    
    return sampler