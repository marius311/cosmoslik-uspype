#!/usr/bin/env python

import sys, os, pycosmc, pkgutil
import argparse
from pycosmc import params

parser = argparse.ArgumentParser(prog='pycosmc.py')
parser.add_argument('params.ini',nargs='?',help='a parameter file to run')
parser.add_argument('--list',action='store_true',default=False,help='list available modules')
parser.add_argument('--doc',nargs=1,metavar='<module>',help='get the documentation for a module')
parser.add_argument('--build',nargs='?',metavar='<modules>',default=False,help='run build script for a module (default: all modules)')
parser.add_argument('-n',nargs=1,metavar='<# of chains>',default=False,help='run multiple chains with MPI')
parser.add_argument('--qsub',action='store_true',default=False,help='submit via qsub')


if not sys.argv[1:]: parser.print_help()
else:
    
    args = vars(parser.parse_args())
    
    if args['list']:
        pkgs = ['pycosmc.likelihoods', 'pycosmc.derivers', 'pycosmc.models', 'pycosmc.samplers']
        print "Found the following modules:"
        for p in pkgs:
            for _, modname, _ in pkgutil.iter_modules(__import__(p,fromlist=[p.split('.')[1]]).__path__):
                print '  %s.%s'%(p.split('.')[1],modname)
        print "See 'pycosmc.py --doc <module>' for more information on a given module."
        
    elif args['doc']:
        modname = args['doc'][0]
        try:
            mod = __import__('pycosmc.%s'%modname,fromlist=[modname.split('.')[-1]])
        except ImportError:
            print "'%s' module not found.\nSee 'pycosmc.py --list' to list all available modules."%modname
        else:
            print "Documentation for module '%s':"%modname
            print mod.__getattribute__(modname.split('.')[-1]).__doc__
            
    elif args['build'] is not False:
        pycosmc.build(args['build'])
        
    elif args['qsub']:
        from subprocess import Popen, PIPE

        inifile = args['params.ini']
        
        if ('camb' in params.read_ini(inifile).get('models',[])): nnodes, ppn, wall = 6, 1, 24
        else: nnodes, ppn, wall = 1, 8, 6
        
        name = inifile.replace('pycosmc.','').replace('params.','')
        dir = os.path.dirname(os.path.abspath(inifile))
        sys.argv.remove('--qsub')
        
        proc = Popen(["qsub"," -q","usplanck","-l",
                      "nodes=%s:ppn=%s,pvmem=20gb,walltime=%s:00:00"%(nnodes,ppn,wall),
                      "-N",name,"-o","%s.log"%name,"-j","oe","-V"],stdin=PIPE,stdout=PIPE)
        
        proc.stdin.write('cd %s && %s -m pycosmc %s'%(dir,sys.executable,' '.join(sys.argv[1:])))
        print proc.stdout.readline()

    elif args['n']:
        i = sys.argv.index('-n')
        sys.argv.pop(i); sys.argv.pop(i)
        os.system("mpiexec -n %i %s -m pycosmc %s"%(int(args['n'][0]),sys.executable,' '.join(sys.argv[1:])))
        
    elif args['params.ini']:
        for _ in pycosmc.pycosmc(args['params.ini']): pass

