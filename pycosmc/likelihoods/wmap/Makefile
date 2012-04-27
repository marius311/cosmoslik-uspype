include ../../../Makefile.inc

pywmap:  pywmap.f90
	f2py --f90flags=-I$(WMAPDIR) -c pywmap.f90 $(LAPACK) -L$(WMAPDIR) -lwmap7 $(CFITSIO) -m pywmap $(F2PYFLAGS)

clean:
	-rm -f pywmap.so
