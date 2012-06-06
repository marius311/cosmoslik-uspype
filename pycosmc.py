#!/usr/bin/env python

import sys, os, pycosmc, pkgutil
import argparse

if __name__=="__main__":
    
    parser = argparse.ArgumentParser(prog='pycosmc.py')
    parser.add_argument('params.ini',nargs='?',help='a parameter file to run')
    parser.add_argument('--list',action='store_true',default=False,help='list available modules')
    parser.add_argument('--doc',nargs=1,metavar='<module>',help='get the documentation for a module')
    
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
                mod = __import__('pycosmc.%s'%modname,fromlist=[modname.split('.')[1]])
            except ImportError:
                print "'%s' module not found. See 'pycosmc.py --help' to list all available modules."%modname
            else:
                print "Documentation for module '%s':"%modname
                print mod.__getattribute__(modname.split('.')[1]).__doc__
        elif args['params.ini']:
            for _ in pycosmc.pycosmc(args['params.ini']): pass
