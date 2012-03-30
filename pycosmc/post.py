import os, sys, re
from numpy import *
from matplotlib.pyplot import *
from matplotlib.mlab import movavg

class Chain(dict):
    """
    An MCMC chain. This is just a dictionary mapping parameter names
    to lists of values, along with the special keys 'lnl' and 'weight'
    """
    def __init__(self,*args,**kwargs):
        super(Chain,self).__init__(*args,**kwargs)
        if 'weight' not in self: self['weight']=ones(len(self.values()[0]))
        
    def params(self): 
        """Returns the parameters in this chain (i.e. the keys except 'lnl' and 'weight'"""
        return set(self.keys())-set(["lnl","weight"])
    
    def sample(self,s,keys=None): 
        """Return a sample or a range of samples depending on if s is an integer or a slice object."""
        return Chain((k,self[k][s]) for k in (keys if keys else self.keys()))
    
    def matrix(self,params=None):
        """Return this chain as an nsamp * nparams matrix."""
        return vstack([self[p] for p in (params if params else self.params())]).T
    
    def cov(self,params=None): 
        """Returns the covariance of the parameters (or some subset of them) in this chain."""
        return get_covariance(self.matrix(params), self["weight"])
    
    def mean(self,params=None): 
        """Returns the mean of the parameters (or some subset of them) in this chain."""
        return average(self.matrix(params),axis=0,weights=self["weight"])
    
    def std(self,params=None): 
        """Returns the std of the parameters (or some subset of them) in this chain."""
        return sqrt(average((self.matrix(params)-self.mean(params))**2,axis=0,weights=self["weight"]))
    
    def acceptance(self): 
        """Returns the acceptance ratio."""
        return 1./mean(self["weight"])
    
    def thin(self,delta):
        """Take every delta samples."""
        c=ceil(cumsum([0]+self['weight'])/float(delta))
        ids=where(c[1:]>c[:-1])[0]
        weight=diff(c[[0]+list(ids+1)])
        t=self.sample(ids)
        t['weight']=weight
        return t
    
    def savecov(self,file,params=None):
        """Write the covariance to a file where the first line is specifies the parameter names."""
        if not params: params = self.params()
        with open(file,'w') as f:
            f.write("# "+" ".join(params)+"\n")
            savetxt(f,self.cov(params))
            
    def savechain(self,file,params=None):
        """Write the chain to a file where the first line is specifies the parameter names."""
        keys = ['lnl','weight']+list(params if params else self.params())
        with open(file,'w') as f:
            f.write("# "+" ".join(keys)+"\n")
            savetxt(f,self.matrix(keys))
            
    def plot(self,param):
        """Plot the value of a parameter as a function of sample number."""
        plot(cumsum(self['weight']),self[param])
        
    def like1d(self,p,**kw): 
        """Plots 1D likelihood contours for a parameter."""
        likelihoodplot1d(self[p],weights=self["weight"],**kw)
        
    def like2d(self,p1,p2,**kw): 
        """Plots 2D likelihood contours for a pair of parameters."""
        likelihoodplot2d(self[p1], self[p2], weights=self["weight"], **kw)

        
        
class Chains(list):
    """A list of chains, probably from several MPI runs"""
    
    def burnin(self,nsamp): 
        """Remove the first nsamp samples from each chain."""
        return Chains(c.sample(slice(nsamp,-1)) for c in self)
    
    def join(self): 
        """Combine the chains into one."""
        return Chain((k,hstack([c[k] for c in self])) for k in self[0].keys())
    
    def plot(self,param): 
        """Plot the value of a parameter as a function of sample number for each chain."""
        for c in self: c.plot(param)
    
def likelihoodplot2d(datx,daty,weights=None,nbins=15,which=[.68,.95],filled=False,color='k',**kw):
    if (weights==None): weights=ones(len(datx))
    H,xe,ye = histogram2d(datx,daty,nbins,weights=weights)
    xem, yem = movavg(xe,2), movavg(ye,2)
    (contourf if filled else contour)(xem,yem,transpose(H),levels=confint2d(H, which[::-1]+[0]),colors=color,**kw)
    
def likelihoodplot1d(dat,weights=None,nbins=30,range=None,maxed=True,**kw):
    if (weights==None): weights=ones(len(dat))
    H, xe = histogram(dat,bins=nbins,weights=weights,normed=True,range=range)
    if maxed: H=H/max(H)
    xem=movavg(xe,2)
    plot(xem,H,**kw)

def get_covariance(data,weights=None):
    if (weights==None): return cov(data.T)
    else:
        mean = sum(data.T*weights,axis=1)/sum(weights)
        zdata = data-mean
        return dot(zdata.T*weights,zdata)/(sum(weights)-1)

def confint2d(hist,which):
    """Return """
    H=sort(hist.ravel())[::-1]
    sumH=sum(H)
    cdf=array([sum(H[H>x])/sumH for x in H])
    return interp(which,cdf,H)


def load_chain(filename):
    """
    If filename is a chain, return a Chain object.
    If filename is a prefix such that there exists filename_1, filename_2, etc... returns a Chains object
    """
    def load_one_chain(filename):
        with open(filename) as file:
            names = re.sub("#","",file.readline()).split()
            try: data = loadtxt(file)
            except: data = None
            
        return Chain([(name,data[:,i] if data!=None else array([])) for (i,name) in enumerate(names)])
    
    dir = os.path.dirname(filename)
    files = [os.path.join(dir,f) for f in os.listdir('.' if dir=='' else dir) if f.startswith(os.path.basename(filename)+'_') or f==os.path.basename(filename)]
    if len(files)==1: return load_one_chain(files[0])
    elif len(files)>1: return Chains(filter(lambda c: c!={}, (load_one_chain(f) for f in files)))
    else: raise IOError("File not found: "+filename) 


