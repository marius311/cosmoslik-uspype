#!/usr/bin/env python

import sys, os, pycosmc

if __name__=="__main__":
    
    if (len(sys.argv) != 2) or sys.argv[1]=='--help':
        print "Usage: python pycosmc.py parameter_file.ini"
    else:
        for _ in pycosmc.pycosmc(sys.argv[1]): pass
