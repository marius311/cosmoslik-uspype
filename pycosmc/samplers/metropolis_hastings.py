from numpy import log, mean, array, sqrt, diag, genfromtxt, sum, dot, cov
from random import random
from numpy.random import multivariate_normal
import pycosmc.mpi as mpi, re, time
from itertools import product
from collections import namedtuple
from pycosmc.modules import Sampler
    
sampletuple = namedtuple('sampletuple',['x','weight','lnl','extra'])
    

class metropolis_hastings(Sampler):
    
    def init(self,p):
        if mpi.get_rank()>0: 
            if 'output_file' in p: p['output_file']+=('_%i'%mpi.get_rank())
            if 'samples' in p: p['samples']/=(mpi.get_size()-1)
            p['quiet']=True 
        elif mpi.get_size()>1: 
            p.pop('output_file')
        
        
    def sample(self,x,lnl,p):
        """
        Returns a generator which yields samples from the likelihood 
        using the Metropolis-Hastings algorithm.
        
        The samples returned are tuples of (x,weight,lnl,extra)
            lnl - likelihood
            weight - the statistical weight (could be 0 for rejected steps)
            x - the vector of parameter values
            extra - the extra information returned by lnl
        """
        
        if mpi.get_size()==1: return _mcmc(x,lnl,p)
        else: return _mpi_mcmc(x,lnl,p)
    
    
    
def _mcmc(x,lnl,p):
    
    (cur_lnl, cur_extra), cur_x, cur_weight = lnl(x,p), x, 1
    
    def update(v): 
        if v!=None: p.update(v)
        
    for _ in range(p.get("samples",100)):
        test_x = multivariate_normal(cur_x,p['_cov'])
        test_lnl, test_extra = lnl(test_x, p)
                
        if (log(random()) < cur_lnl-test_lnl):
            update((yield(sampletuple(cur_x, cur_weight, cur_lnl, cur_extra))))
            cur_lnl, cur_weight, cur_x, cur_extra = test_lnl, 1, test_x, test_extra
        else:
            if p.get('metropolis_hastings',{}).get('weighted_samples',True): cur_weight+=1
            else: update((yield(sampletuple(cur_x, cur_weight, cur_lnl, cur_extra))))
            
            if p.get('metropolis_hastings',{}).get('rejected_samples',True): 
                update((yield(sampletuple(test_x, 0, test_lnl, test_extra))))
            
            
            
def _mpi_mcmc(x,lnl,p):  

    (rank,size,comm) = mpi.get_mpi()
    from mpi4py import MPI #FIX THIS
    
    if rank==0:
        finished = [False]*(size-1)
        samples, weights, lnls = [[[] for _ in range(size-1)] for __ in range(3)] 
        while (not all(finished)):
            while not comm.Iprobe(source=MPI.ANY_SOURCE, tag=0): time.sleep(.1) #Hack so OpenMPI doesn't eat 100% CPU
            (source,new_samples)=comm.recv(source=MPI.ANY_SOURCE)
            if (new_samples!=None): 
                lnls[source-1]+=[s.lnl for s in new_samples]
                samples[source-1]+=[s.x for s in new_samples]
                weights[source-1]+=[s.weight for s in new_samples]
                
                if (p.get("proposal_update",True) 
                    and sum(weights[source-1])>p.get('proposal_update_start',1000)):
                    comm.send({"_cov":get_new_proposal(samples,weights)},dest=source)
                else: comm.send({},dest=source)
                comm.send(None,dest=source)

            else: 
                finished[source-1]=True
                
    else:
        samples = []
        sampler = _mcmc(x, lnl, p)
        s=sampler.next()
        while True:
            try:
                yield s
                samples.append(s)
                if len(samples)==50:
                    comm.send((rank,samples),dest=0)
                    s = sampler.send(comm.recv(source=0))
                    samples = []
                else: s = sampler.next()
            except StopIteration: break
        comm.send((rank,None),dest=0)

       
def get_covariance(data,weights=None):
    if (weights==None): return cov(data.T)
    else:
        mean = sum(data.T*weights,axis=1)/sum(weights)
        zdata = data-mean
        return dot(zdata.T*weights,zdata)/(sum(weights)-1)


def get_new_proposal(samples,weights=None):
    """
    shape(samples) = (nchains,nsamples,nparams)
    shape(weights) = (nchains,nsamples)
    """
    data = array(reduce(lambda a,b: a+b,[s[len(s)/2:] for s in samples]))
    weights = array(reduce(lambda a,b: a+b,[w[len(w)/2:] for w in weights]))
    return get_covariance(data, weights)
    
    
def gelman_rubin_R(samples):
    return 1
    
