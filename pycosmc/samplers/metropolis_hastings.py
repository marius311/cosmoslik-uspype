import sys, os, re
import numpy as np
from numpy import *
from numpy.linalg import inv
from random import random
from itertools import product, repeat
from matplotlib.mlab import movavg
from matplotlib.pyplot import plot, hist, contour, contourf
from ast import literal_eval
from numpy.ma.core import transpose, sort
from numpy.lib.function_base import interp
import mpi 
    
    
def mcmc(start,lnl):
    if mpi.get_size()==1: return _mcmc(start,lnl)
    else: return _mpi_mcmc(start,lnl)


def _mcmc(start,lnl,step_fn=None):
    """
    Run an MCMC chain. 
    
    Parameters:
        start - dictionary of parameter key/values or a filename
        lnl - a function (or list of functions) of (dict:params) which returns the negative log-likelihood
    
    Returns:
        A dictionary of samples. There is a key for every parameter which was varied, 
        one for every derived parameter, and one for "lnl" and "weight"
        
    Example:
        samples = mcmc({"x":"0 [-10 10 1]",'samples':10000}, lambda p: p["x"]**2/2)
        
        plot(samples["x"])
    """
    
    initialize_covariance(start)
    start["$COV_FAC"]=float(start.get("proposal_scale",2.4))

    cur_params = start.copy()
    #Initialize starting sample
    samples = {"lnl":[lnl(cur_params)],"weight":[1]}
    for name in get_outputted(cur_params): samples[name]=[cur_params[name]]
    
    #Initialize file if were writing to a file
    if (cur_params.get("chain","")!=""):
        output_file = open(cur_params["chain"],"w")
        output_file.write("# lnl weight "+" ".join(get_outputted(cur_params))+"\n")
        output_to_file = True
    else:
        output_to_file = False
    
    #Start the MCMC
    print "Starting chain..."
    for _ in range(start.get("samples",100)):
        test_params = dict(start,**propose_step_gaussian(cur_params))
        
        # Check min/max bounds
        test_lnl = 0
        for name in get_varied(test_params): 
            if (not (test_params["*"+name][1] < test_params[name] < test_params["*"+name][2])): test_lnl = np.inf
        
        #Get likelihood
        if (test_lnl != np.inf): test_lnl = lnl(test_params)
                
        if not test_params.get('$MPI',False) and test_params.get("mcmc_verbose",True): 
            print "Like=%.2f Ratio=%.3f Sample=%s" % (test_lnl,1./np.mean(array(samples["weight"])),dict([(name,test_params.get(name)) for name in get_outputted(cur_params)])) 

        if (log(random()) < samples["lnl"][-1]-test_lnl):

            #Add to file (which lags samples by one accepted sample so we get the weight right)
            if (output_to_file): 
                output_file.write(" ".join([str(samples[name][-1]) for name in ["lnl","weight"]+get_outputted(cur_params)])+"\n") 
                output_file.flush()
            
            
            #Add to samples
            for name in get_outputted(test_params): samples[name].append(test_params[name])
            samples["lnl"].append(test_lnl)
            samples["weight"].append(1)
            
            cur_params.update({k:test_params[k] for k in get_varied(cur_params)}) 

        else:
            samples["weight"][-1] += 1
            
            
        #Call step function provided by user
        if step_fn: step_fn(cur_params,samples)
            
    if (output_to_file): output_file.close()
    
    for k in samples.keys(): samples[k]=array(samples[k])

    return samples



def _mpi_mcmc(start,lnl):  
    """
    
    Runs an MCMC chain on an MPI cluster using (processes-1) workers and 1 master.
    
    Example:
    
        mpiexample.py:
            import mcmc
            mpi_mcmc({"*x":[0,-10,10,1],'samples':10000,"file_root":"chain"}, lambda p: p["x"]**2/2)
    
        mpiexec -n 8 python mpiexample.py
    
    
    Internally, every delta_send_samples number of steps, the worker processes communicate new steps to the master process. 
    They then immediately ask the master process for a dictionary of parameters with which to update the local copy (for example, a new proposal covariance).
    
    The master process aggregates the samples from all of the chains, calculating a proposal covariance and sending it out to the
    worker processes whenever necessary. 
    
    """

    from mpi4py import MPI
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()-1

    start["samples"]/=size
    if "chain" in start: start["chain"]+=("_"+str(rank))
    start["$MPI"]=True
    start["$UPDATECOV"]=start.get("proposal_update",True)
    [start["$RMIN"],start["$RMAX"]] = start.get("proposal_update_minmax_R",[0,np.inf])

    delta_send_samples = 2

    def wslice(samples,delta):
        """
        Gets the last delta non-unique samples, possibly slicing one weighted sample
        """
        weights = samples["weight"]
        nws = len(filter(lambda x: x,np.cumsum(weights[-delta:])-sum(weights[-delta:])>-delta))
        sl = dict([(name,samples[name][-nws:]) for name in samples.keys()])
        sl["weight"][0]-=(sum(weights[-nws:])-delta)
        return sl
    
        
    def mpi_step_fn(params,samples):
        """ 
        This is the worker process code.
        """
        if (sum(samples["weight"])%delta_send_samples==0):
            print "Chain %i: steps=%i approval=%.3f best=%.2f" % (rank,sum(samples["weight"]),1./np.mean(array(samples["weight"])),min([inf]+samples["lnl"][1:]))
            comm.send((rank,wslice(samples,delta_send_samples)))  
            new_params = comm.recv()
            if (new_params!=None):  
                print "Chain %i: New proposal: %s"%(rank,zip(get_varied(params),sqrt(diag(new_params["$COV"]))))
                params.update(new_params)
            

    def mpi_master_fn(samples,source):
        """ 
        
        This is the master process code.
        
        samples - a list of samples aggregated from all of the chains
        source - the rank of the chain which is requesting updated parameters
        
        """
        
        if (start["$UPDATECOV"] and (lambda x: x>=2000 and x%delta_send_samples==0)(sum(samples[source-1]["weight"])) and start["$RMIN"]<gelman_rubin_R(samples)<start["$RMAX"]): 
            return {"$COV":get_new_proposal(samples, start)}
        else: 
            return None


    if (rank==0):
        samples = [dict([(name,[]) for name in get_varied(start)+["lnl","weight"]]) for _ in range(size)]
        finished = [False]*size
        while (not all(finished)):
            (source,new_samples)=comm.recv(source=MPI.ANY_SOURCE)
            if (new_samples!=None): 
                for name in get_varied(start)+['lnl','weight']: 
                    samples[source-1][name]+=(lambda x: x if type(x)==list else [x])(new_samples[name])
            else:
                finished[source-1]=True
                
            comm.send(mpi_master_fn(samples,source),dest=source)
        
    else:
        
        _mcmc(start,lnl,mpi_step_fn)
        comm.send((rank,None))


def get_varied(params): return params["$VARIED"]
def get_outputted(params): return params["$OUTPUT"]

def propose_step_gaussian(params,fac=None):
    """Take a gaussian step in $VARIED according to $COV"""
    varied_params = get_varied(params)
    cov = params["$COV"]
    nparams = len(varied_params)
    if (shape(cov)!=(nparams,nparams)):
        raise ValueError("Covariance not the same length as number varied parameters.")
    dxs = np.random.multivariate_normal([0]*nparams,cov) * sqrt(fac if fac else params["$COV_FAC"])
    propose = params.copy()
    return {name:params[name]+dx for (name,dx) in zip(varied_params,dxs)}

def initialize_covariance(params):
    """
    Load the covariance, defaulting to diagonal entries from the WIDTH of each parameter
    """
    
    v=params.get("proposal_matrix","")
    if (v==""): prop_names, prop = [], None
    else: 
        with open(v) as f:
            prop_names = re.sub("#","",f.readline()).split()
            prop = genfromtxt(f)

    params["$COV"] = np.diag([params["*"+name][3]**2 for name in get_varied(params)])
    common = set(get_varied(params)) & set(prop_names)
    if common: 
        idxs = zip(*(list(product([ps.index(n) for n in common],repeat=2)) for ps in [get_varied(params),prop_names]))
        for ((i,j),(k,l)) in idxs: params["$COV"][i,j] = prop[k,l]
   
def get_covariance(data,weights=None):
    if (weights==None): return np.cov(data.T)
    else:
        mean = np.sum(data.T*weights,axis=1)/np.sum(weights)
        zdata = data-mean
        return np.dot(zdata.T*weights,zdata)/(np.sum(weights)-1)

def get_new_proposal(samples,params):
    nchains = len(samples)
    nsamp = [len(s['weight']) for s in samples]
    data = array([reduce(lambda a,b: a+b,[samples[i][name][:nsamp[i]/2] for i in range(nchains)]) for name in get_varied(params)]).T
    weights = array(reduce(lambda a,b: a+b,[samples[i]["weight"][:nsamp[i]/2] for i in range(nchains)]))
    return get_covariance(data, weights)
    
def gelman_rubin_R(samples):
    return 1
    

