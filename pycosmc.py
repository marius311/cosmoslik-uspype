#!/usr/bin/env python

import sys, os, pycosmc

if __name__=="__main__":
    
    if (len(sys.argv) != 2) or sys.argv[1]=='--help':
        print "Usage: python pycosmc.py parameter_file.ini"
        for m in ['likelihoods','models','derivers','samplers']:
            print 'Available %s:'%m
            for f in os.listdir(os.path.join(os.path.dirname(__file__),'pycosmc',m)):
                try: print ' %s'%__import__('pycosmc.%s.%s'%(m,f),fromlist=f).__name__
                except: pass
    else:
        for _ in pycosmc.pycosmc(sys.argv[1]): pass
