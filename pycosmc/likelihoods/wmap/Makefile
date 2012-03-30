include ../../Makefile.inc

pywmap:  pywmap.f90
	f2py --f90flags=-I$(WMAPDIR) -c pywmap.f90 --link-lapack_opt -L$(WMAPDIR) -lwmap7 $(CFITSIO) -m pywmap $(F2PYFLAGS)
