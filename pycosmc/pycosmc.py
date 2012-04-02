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
    loadmods = lambda x1: [__import__('pycosmc.%s.%s'%(x1,m),fromlist=m) for m in p.get(x1,'').split()]
    lnls, models, derivers, samplers = [loadmods(x1) for x1 in ['likelihoods','models','derivers','samplers']]
    
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
    samples = namedtuple('sampletuple',['x','weight','lnl','params'])([],[],[],[])
    for sampler in samplers:
        print "Starting %s sampler..."%sampler.__name__.split('.')[-1]
        for (nsamp,s) in enumerate(sampler.sample([p[k] for k in p['$SAMPLED']],lnl,**p)):
            yield s
            x1, w1, l1, p1 = s
                          
            #Add derived if they're not in there
            if p1==None or not all(k in p1 for k in p['$OUTPUT']): 
                p1 = p.copy(); p.update(p1); p.update(zip(p['$SAMPLED'],x1))
                for d in derivers: d.add_derived(p1)

            if w1!=0: 
                for (l,v) in zip(samples,(x1, w1, l1, p1)): l.append(v)

            if w1!=0 and f!=None: 
                f.write(' '.join(map(str,[l1,w1]+[p1.get(name,'nan') for name in p['$OUTPUT']]))+'\n')
                f.flush()
                
            if nsamp%p.get('update_frequency',1)==0:
                print "%ssamples=%s best=%.3g acceptance=%.3f mean={%s} proposal={%s}" % \
                    ('' if mpi.get_rank()==0 else 'Chain %i: '%mpi.get_rank(),
                     str(sum(samples.weight)),
                     min(samples.lnl+[inf]),
                     1./mean(samples.weight),
                     ', '.join(['%s:%.3g'%(k,v) for k,v in dict([(name,mean([e[name] for e in samples.params[len(samples.params)/2:]])) for name in p['$OUTPUT']]).items()]),
                     (', '.join(['%s:%.3g'%(k,v) for k,v in dict(zip(p['$SAMPLED'],sqrt(diag(p['$COV'])))).items()])) if '$COV' in p else ''
                     ) 
        
    if f!=None: f.close()