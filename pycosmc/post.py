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
        
    def copy(self):
        """Deep copy the chain so post-processing, etc... works right"""
        return Chain({k:v.copy() for k,v in self.iteritems()})
        
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
        like1d(self[p],weights=self["weight"],**kw)
        
    def like2d(self,p1,p2,**kw): 
        """Plots 2D likelihood contours for a pair of parameters."""
        like2d(self[p1], self[p2], weights=self["weight"], **kw)
        
    def likegrid(self,**kwargs):
        """Plot several 2d likelihood contours."""
        likegrid(self, **kwargs)

        
        
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
    
def like2d(datx,daty,weights=None,nbins=15,which=[.68,.95],filled=False,color='k',**kw):
    if (weights==None): weights=ones(len(datx))
    H,xe,ye = histogram2d(datx,daty,nbins,weights=weights)
    xem, yem = movavg(xe,2), movavg(ye,2)
    (contourf if filled else contour)(xem,yem,transpose(H),levels=confint2d(H, which[::-1]+[0]),colors=color,**kw)
    
def like1d(dat,weights=None,nbins=30,range=None,maxed=True,**kw):
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


def likegrid(chains,ps=None,fig=None,colors=['b','g','r'],nbins1d=30, nbins2d=20):
    if fig==None: fig=figure()
    if type(chains)!=list: chains=[chains]
    if ps==None: ps = sorted(reduce(lambda x,y: set(x)&set(y), [c.params() for c in chains]))
    colors=colors[:len(chains)]
    #figsize(20,20)
    subplots_adjust(hspace=0,wspace=0)

    c=chains[0]
    lims = [(max(min(c[p]),mean(c[p])-4*std(c[p])),min(max(c[p]),mean(c[p])+4*std(c[p]))) for p in ps]
    ticks = [[t for t in ts if l[0]<=t<=l[1]] for (ts,l) in zip((c.mean(ps)+c.std(ps)*transpose([[-2,0,2]])).T ,lims)]
    n=len(ps)
    axs=n*[n*[None]]
    for (i,p1) in enumerate(ps):
        for (j,p2) in enumerate(ps):
            if (i<=j):
                ax=axs[i][j]=subplot(n,n,j*n+i+1)
                xlim(*lims[i])
                ax.set_xticks(ticks[i])
                if (i==j): 
                    for (ch,col) in zip(chains,colors): 
                        if p1 in ch: ch.like1d(p1,nbins=nbins1d,color=col)
                    ax.set_yticks([])
                    
                elif (i<j): 
                    for (ch,col) in zip(chains,colors): 
                        if p1 in ch and p2 in ch: ch.like2d(p1,p2,filled=False,nbins=nbins2d,color=col)
                    ylim(*lims[j])
                    ax.set_yticks(ticks[j])
                        
                if i==0: 
                    ylabel(p2)
                    ax.set_yticklabels(['%.3g'%t for t in ticks[j]])
                else: 
                    ax.set_yticklabels([])
                
                if j==n-1: 
                    xlabel(p1)
                    ax.set_xticklabels(['%.3g'%t for t in ticks[i]])
                else: 
                    ax.set_xticklabels([])
    return axs


def confint2d(hist,which):
    """Return """
    H=sort(hist.ravel())[::-1]
    sumH=sum(H)
    cdf=array([sum(H[H>x])/sumH for x in H])
    return interp(which,cdf,H)


def load_chain(path,paramnames=None):
    """
    If path is a chain, return a Chain object.
    If path is a prefix such that there exists path_1, path_2, etc... returns a Chains object
    """
    def load_one_chain(path):
        if os.path.isdir(path):
            if paramnames!=None: raise Exception("Can't specify custom parameter names if loading chain from a directory.")
            chain = {}
            for k in os.listdir(path):
                try: chain[k]=loadtxt(os.path.join(path,k),usecols=[-1])
                except: pass
            return Chain(chain)
        else:
            names = None
            if paramnames==None:
                pnfiles = [os.path.join(os.path.dirname(path),f) for f in os.listdir(os.path.dirname(path)) if f.endswith('.paramnames') and os.path.basename(path).startswith(f[:-len('.paramnames')])]
                if len(pnfiles)>1: raise Exception('Found multiple paramnames files for this chain; %s'%pnfiles)
            
            if paramnames or pnfiles:
                with open(paramnames or pnfiles[0]) as f:
                    names = ['weight','lnl']+[line.split()[0] for line in f]
                    
            with open(path) as f:
                if names==None: names = re.sub("#","",f.readline()).split()
                try: data = loadtxt(f).T
                except: data = [array([])]*len(names)
                
            return Chain(zip(names,data))
    
    path=os.path.abspath(path)
    dir = os.path.dirname(path)
    files = [os.path.join(dir,f) for f in os.listdir('.' if dir=='' else dir) if re.match(os.path.basename(path)+'_[0-9]+',f) or f==os.path.basename(path)]
    if len(files)==1: return load_one_chain(files[0])
    elif len(files)>1: return Chains(filter(lambda c: c!={}, (load_one_chain(f) for f in files)))
    else: raise IOError("File not found: "+path) 


