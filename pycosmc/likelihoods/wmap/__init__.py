"""
===============
WMAP Likelihood
===============

- Written by WMAP team (see `<http://lambda.gsfc.nasa.gov/>`_)
- CosmoSlik module by Marius Millea
- Updated July 1, 2012

Description
===========
    
This module wraps the official WMAP likelihood code. 
Some minor modifications were made to allow:

- Choosing the WMAP data directory at runtime
- Choosing the lmin/lmax at runtime


Install Notes
=============

This build this module run::

    ./pycosmc.py --build likelihoods.wmap
    
The Makefile for this module reads the following flags from ``Makefile.inc``:

- ``$(CFITSIO)``
- ``$(LAPACK)``
- ``$(F2PYFLAGS)``


Models
======

The WMAP module requires a `Model` which provides the following:

- ``cl_TT``
- ``cl_TE``
- ``cl_EE``
- ``cl_BB``
    
Extra-galactic foregrounds are ignored. 


Parameters
==========

This module reads the following parameters from the ini file:

[wmap].data_dir
---------------
    The path to the wmap/data directory.

[wmap].use
----------
    A subset of ``['TT','TE','EE','BB']`` corresponding to 
    which likelihood terms to use.
    
[wmap].TT.lrange
----------------
    The TT range in ell to use in the likelihood

[wmap].TE.lrange
----------------
    The TE range in ell to use in the likelihood

"""

from wmap import wmap

