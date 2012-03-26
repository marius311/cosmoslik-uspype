#!/usr/bin/env python

import sys
from ini import read_ini
from operator import __add__

lnls = []
models = []
derivers = []

def lnl(p,derivative=0):
    #Calculate derived parameters
    for d in derivers: d.add_derived(p)
    
    #Evaluate models
    model = reduce(lambda x,y: dict(x,**y),(m.get(p) for m in models),{})
    
    #Evaluate likelihoods
    elnls = (l.lnl(model,p,derivative=derivative) for l in lnls)

    return sum(elnls) if derivative==0 else map(__add__, *elnls)


def pycosmc(p):
    if isinstance(p,str): p=read_ini(p)

    #Load model, likelihood, and derivers parameter modules
    global lnls, models, derivers
    loadmod = lambda x: [__import__(x+'.'+l,fromlist=l) for l in p[x].split()]
    lnls, models, derivers = [loadmod(x) for x in ['likelihoods','models','derivers']]
    
    #Initialize modules
    for i in models+lnls+derivers:
        if 'init' in i.__dict__: i.init(p)
        
    #Load and call sampler
    mcmc = __import__('samplers.'+p['sampler'], fromlist=p['sampler'])
    return mcmc.mcmc(p,lnl)
    
    
if __name__=="__main__":
    
    if (len(sys.argv) != 2): 
        print "Usage: python signal_to_params.py parameter_file.ini"
        sys.exit()

    pycosmc(sys.argv[1])

    
    
