from matplotlib.pyplot import *
from numpy import *
import os

def lnl(model,p,**kw):
    ion()
    clf()
    ax=subplot(111)
    for i in p['_likelihoods']+p['_models']+p['_derivers']+p['_samplers']:
        if 'diagnostic' in i.__dict__: i.diagnostic({'cl_TT':ax},p)
    ax.set_yscale('log')
    draw()
    ioff()
    return 0