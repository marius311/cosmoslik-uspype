#!/usr/bin/env python

import sys, pycosmc, pkgutil

if __name__=="__main__":
    
    if (len(sys.argv) == 2): 
        if sys.argv[2]=='help':
            print "Usage: python pycosmc.py parameter_file.ini"
            
            
        
        sys.exit()
        
    for _ in pycosmc.pycosmc(sys.argv[1]): pass
