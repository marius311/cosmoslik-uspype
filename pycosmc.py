#!/usr/bin/env python

import sys, pycosmc

if __name__=="__main__":
    
    if (len(sys.argv) != 2): 
        print "Usage: python pycosmc.py parameter_file.ini"
        sys.exit()
        
    for _ in pycosmc.pycosmc(sys.argv[1]): pass
