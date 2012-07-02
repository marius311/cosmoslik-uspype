"""

====
CAMB
====

- Written by Antony Lewis and Anthony Challinor (see `<http://camb.info>`_)
- CosmoSlik module by Marius Millea


Usage
=====

To use this module, add ``camb`` to the list of ``models``.

Notes
=====

CosmoSlik uses a modified CAMB which under-the-hood reads the ini file 
from `stdin` and writes results to `stdout`. This module also handles some
of the logic in chosing settings. In particular,

- It automatically choses which of ``get_scalar``, ``get_tensor``, 
  ``get_vector``, and ``get_transfer`` should be set to true based 
  on which likelihoods are used. 
- It automatically sets ``do_lensing=True`` if ``Alens>0``.
- If automatically increases ``lmax`` used in calculations to ensure 
  full accuracy of outputs up to the specified ``lmax``. 

.. note:: 

    For users wishing to modify or add functionality to CAMB, we recommend starting
    by making a copy of this module and leaving ``inidriver.F90`` mostly
    intact (it contains the under-the-hood changes mentioned above.)

Parameters
==========

[camb].X
--------
    X can be any valid CAMB parameter that you would normally give in the
    CAMB ini file. These parameters are passed straight to CAMB (and can
    be sampled over).  
    
[camb].defaults
---------------
    Path to a CAMB ini file which specifies default values for parameters not set 
    by the user or by CosmoSlik. (default: defaults.ini in the module folder) 
    
Alens
-----
    cl = Alens * cl_lensed + (1-Alens) * cl_scalar 


"""

from camb import camb