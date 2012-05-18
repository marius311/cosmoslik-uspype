from operator import __add__
from numpy import mean, sqrt, diag, inf, std, loadtxt
from collections import namedtuple
from itertools import product, chain
import mpi, re
import params

derivers = []
models = []
lnls = []
samplers = []

def lnl(x,p):
    
    #Convert vector x to nice named dictionary
    p = p.copy(); p.update(zip(p.get_all_sampled().keys(),x))
    
    #Calculate derived parameters
    for d in derivers: d.add_derived(p)
    
    #Check priors
    if not all(v[1] < p[k] < v[2] for k, v in p.get_all_sampled().items()): 
        return inf, p
    else: 
        #Evaluate models and call likelihoods
        model = ModelDict()
        for m in models: model.update(m.get(p,p['_model.required']))
        return (sum(l.lnl(p,model.for_module(l)) for l in lnls),p)


def pycosmc(p,**kwargs):
    p=params.read_ini(p) if isinstance(p,str) else p
    p.update(kwargs)
    params.eval_values(p)
    params.process_parameters(p)

    #Import the various modules
    for (l,k) in zip((lnls,models,derivers,samplers),('likelihoods','models','derivers','samplers')):
        x = p.get(k,[])
        if isinstance(x,str): x=x.split()
        l += [__import__('pycosmc.%s.%s'%(k,m),fromlist=m).__getattribute__(m)() if isinstance(m,str) else m for m in x]  
     

    #Initialize modules
    for i in lnls+models+derivers+samplers:
        if hasattr(i,'init'): 
            print 'Initializing %s...'%i.__class__.__name__
            i.init(p)
            
    p['_model.required']=set(chain(*[l.get_required_models(p) for l in lnls]))
    
    #Sampled and outputted parameters and covariance   
    for l in lnls:
        for k, v in l.get_extra_params(p).items(): 
            p.setdefault(l.__class__.__name__).add_sampled_param(k,*v)
    sampled = p.get_all_sampled().keys()
    outputted = sampled + [(k,) for k in p.get('derived','').split()]
    p['_cov'] = initialize_covariance(p)
    
    #Prep output file
    if 'output_file' in p: 
        f = open(p['output_file'],'w')
        f.write("# lnl weight "+" ".join(['.'.join(k) for k in outputted])+"\n")
    else: f = None
        
    #Run samplers
    samples = namedtuple('sampletuple',['x','weight','lnl','params'])([],[],[],[])
    for sampler in samplers:
        print "Starting %s sampler..."%sampler.__class__.__name__
        for (nsamp,s) in enumerate(sampler.sample([p[k] for k in sampled],lnl,p),1):
            yield s
            x1, w1, l1, p1 = s
                          
#            #Add derived if they're not in there
#            if p1==None or not all(k in p1 for k in p['_output']): 
#                p1 = p.copy(); p.update(p1); p.update(zip(p['_sampled'],x1))
#                for d in derivers: d.add_derived(p1)
#                assert all(k in p1 for k in outputted), "Derivers didn't calculate all the derived parameters. Check 'output' key or add derivers."

            if w1!=0:
                for (l,v) in zip(samples,(x1, w1, l1, p1)): l.append(v)

            if f!=None and w1!=0: 
                f.write(' '.join(map(str,[l1,w1]+[p1[name] for name in outputted]))+'\n')
                f.flush()
                
            if nsamp%p.get('update_frequency',1)==0:
                print "%saccepted=%s/%i(%.1f%%) best=%.2f last={%s}" % \
                    ('' if mpi.get_rank()==0 else 'Chain %i: '%mpi.get_rank(),
                     len(samples.weight),
                     nsamp,
                     100*float(len(samples.weight))/nsamp,
                     min(samples.lnl+[inf]),
                     ', '.join([('like:%.2f'%l1)]+['%s:%.4g'%('.'.join(name),p1.get(name,float('nan'))) for name in outputted])
                     ) 

    if f!=None: f.close()
    
    if 'dump_samples' in p:
        import cPickle
        with open(p['dump_samples'],'w') as f: cPickle.dump(tuple(samples), f, 2)
        
        
class ModelDict(dict):
    
    def __getitem__(self,k):
        try: return dict.__getitem__(self,k)
        except: raise Exception("The likelihood module '%s' needs a model for '%s' but no models provided it. Try running 'pycosmc.py --help' to list available models."%(self._for_module,k))
        
    def for_module(self,l):
        self._for_module = l.__class__.__name__
        return self
        

def initialize_covariance(params):
    """Load the sigma, defaulting to diagonal entries from the WIDTH of each parameter."""
    v=params.get("proposal_matrix","")
    if (v==""): prop_names, prop = [], None
    else: 
        with open(v) as f:
            prop_names = [tuple(k.split('.')) for k in re.sub("#","",f.readline()).split()]
            prop = loadtxt(f)
    sampled = params.get_all_sampled()
    sigma = diag([v[3]**2 for v in sampled.values()])
    common = set(sampled.keys()) & set(prop_names)
    if common: 
        idxs = zip(*(list(product([ps.index(n) for n in common],repeat=2)) for ps in [sampled.keys(),prop_names]))
        for ((i,j),(k,l)) in idxs: sigma[i,j] = prop[k,l]
    return sigma/len(sampled)*params.get('proposal_scale',2.4)**2
