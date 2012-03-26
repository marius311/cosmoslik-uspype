pywmap:  pywmap.f90
	f2py --f90flags=-I$(WMAP) -c pywmap.f90 --link-lapack_opt -I$(WMAP) -L$(WMAP) -lwmap7 -lcfitsio -m pywmap
