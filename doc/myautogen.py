import pkgutil
from textwrap import dedent

pkgs = ['pycosmc.likelihoods', 'pycosmc.derivers', 'pycosmc.models', 'pycosmc.samplers']
for p in pkgs:
    for _, modname, _ in pkgutil.iter_modules(__import__(p,fromlist=[p.split('.')[1]]).__path__):
        print '  %s.%s'%(p.split('.')[1],modname)
        try:
            mod = __import__('%s.%s'%(p,modname),fromlist=[modname.split('.')[-1]])
        except ImportError:
            print "'%s' module not found."%modname
        else:
            doc = mod.__doc__
            if doc!=None: 
                with open('%s.%s'%(p,modname)+'.rst','w') as f: f.write(dedent(doc))


