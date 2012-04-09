from ini import read_ini
from operator import __add__
from numpy import mean, sqrt, diag, inf, std
from collections import namedtuple
import mpi

derivers = []
models = []
lnls = []
samplers = []

def lnl(x,derivative=0,**kwargs):
    
    #Convert vector x to nice named dictionary
    p = kwargs.copy(); p.update(zip(kwargs['_sampled'],x))
    
    #Calculate derived parameters
    for d in derivers: d.add_derived(p)
    
    #Check priors
    if not all(p['*'+k][1] < p[k] < p['*'+k][2] for k in p['_sampled']): elnls = [inf]
    else: 
        #Evaluate models and call likelihoods
        p['_model'] = {}
        for m in models: p['_model'].update(m.get(p))
        elnls = (l.lnl(p['_model'],p,derivative=derivative) for l in lnls)

    return (sum(elnls) if derivative==0 else map(__add__, *elnls)), p


def pycosmc(p,**kwargs):
    p=read_ini(p) if isinstance(p,str) else p
    p.update(kwargs)

    #Import the various modules
    for (l,x) in zip((lnls,models,derivers,samplers),('likelihoods','models','derivers','samplers')): 
        l += [__import__('pycosmc.%s.%s'%(x,m),fromlist=m) for m in p.get(x,'').split()] 
     
    #Likelihood modules will add to this set the quantities they need calculated
    p['_models.get']=set()

    #Initialize modules
    for i in lnls+models+derivers+samplers:
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
        for (nsamp,s) in enumerate(sampler.sample([p[k] for k in p['_sampled']],lnl,**p)):
            yield s
            x1, w1, l1, p1 = s
                          
            #Add derived if they're not in there
            if p1==None or not all(k in p1 for k in p['$OUTPUT']): 
                p1 = p.copy(); p.update(p1); p.update(zip(p['_sampled'],x1))
                for d in derivers: d.add_derived(p1)
                assert all(k in p1 for k in p['$OUTPUT']), "Derivers didn't calculate all the derived parameters. Check 'output' key or add derivers."

            if w1!=0: 
                for (l,v) in zip(samples,(x1, w1, l1, p1)): l.append(v)

            if f!=None: 
                f.write(' '.join(map(str,[l1,w1]+[p1[name] for name in p['$OUTPUT']]))+'\n')
                f.flush()
                
            if nsamp%p.get('update_frequency',1)==0:
                print "%saccepted=%s/%i(%.1f%%) best=%.2f last={%s}" % \
                    ('' if mpi.get_rank()==0 else 'Chain %i: '%mpi.get_rank(),
                     len(samples.weight),
                     nsamp,
                     100./mean(samples.weight),
                     min(samples.lnl+[inf]),
                     ', '.join([('like:%.2f'%l1)]+['%s:%.4g'%(name,p1[name]) for name in p['$OUTPUT']])
                     ) 

    if f!=None: f.close()
    
    if 'dump_samples' in p:
        import cPickle
        with open(p['dump_samples'],'w') as f: cPickle.dump(tuple(samples), f, 2)
