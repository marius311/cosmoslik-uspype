from ini import read_ini
from operator import __add__
from numpy import mean, sqrt, diag, inf
from collections import namedtuple
import mpi

lnls = []
models = []
derivers = []
samplers = []

def lnl(x,derivative=0,**kwargs):
    #Convert vector x to nice named dictionary
    p = dict(kwargs,**dict(zip(kwargs['$SAMPLED'],x)))
    
    #Calculate derived parameters
    for d in derivers: d.add_derived(p)
    
    #Evaluate models
    model = reduce(lambda x,y: dict(x,**y),(m.get(p) for m in models),{})
    
    #Evaluate likelihoods
    elnls = (l.lnl(model,p,derivative=derivative) for l in lnls)

    return (sum(elnls) if derivative==0 else map(__add__, *elnls)), p


def pycosmc(p,**kwargs):
    p=dict(read_ini(p) if isinstance(p,str) else p,**kwargs)

    #Load model, likelihood, and derivers parameter modules
    global lnls, models, derivers, samplers
    loadmods = lambda x: [__import__('pycosmc.%s.%s'%(x,m),fromlist=m) for m in p.get(x,'').split()]
    lnls, models, derivers, samplers = [loadmods(x) for x in ['likelihoods','models','derivers','samplers']]
    
    #Likelihood modules will add to this set the quantities they need calculated
    p['models.calculate']=set()

    #Initialize modules
    for i in models+derivers+lnls+samplers:
        if 'init' in i.__dict__: i.init(p)
        
    #Prep output file
    if 'output_file' in p: 
        f = open(p['output_file'],'w')
        f.write("# lnl weight "+" ".join(p['$OUTPUT'])+"\n")
    else: f = None
        
    #Run samplers
    samples = namedtuple('sampletuple',['x','weight','lnl','extra'])([],[],[],[])
    for sampler in samplers:
        print "Starting %s sampler..."%sampler.__name__.split('.')[-1]
        for (nsamp,s) in enumerate(sampler.sample([p[k] for k in p['$SAMPLED']],lnl,**p)):
            yield s
            x, weight, like, p = s
            if weight!=0: 
                for (l,v) in zip(samples,s): l.append(v)
                            
            if weight!=0 and f!=None: 
                f.write(' '.join(map(str,[like,weight]+[p[name] for name in p['$OUTPUT']]))+'\n')
                f.flush()
                
            if nsamp%p.get('update_frequency',1)==0:
                print "%ssamples=%i best=%.2g acceptance=%.3f mean={%s} proposal={%s}" % \
                    ('' if mpi.get_rank()==0 else 'Chain %i: '%mpi.get_rank(),
                     sum(samples.weight),
                     min(samples.lnl+[inf]),
                     1./mean(samples.weight),
                     ', '.join(['%s:%.3g'%(k,v) for k,v in dict([(name,mean([e[name] for e in samples.extra[len(samples.extra)/2:]])) for name in p['$OUTPUT']]).items()]),
                     ', '.join(['%s:%.3g'%(k,v) for k,v in dict(zip(p['$SAMPLED'],sqrt(diag(p['$COV'])))).items()])
                     ) 
        
    if f!=None: f.close()