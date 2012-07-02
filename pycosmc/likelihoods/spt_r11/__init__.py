"""

========================
SPT Reichardt et al 2011
========================

- Paper: `<http://adsabs.harvard.edu/abs/2011arXiv1111.0932R>`_
- Data: `<http://pole.uchicago.edu/public/data/reichardt11/index.html>`_
- CosmoSlik module by Marius Millea

Usage
=====

To use this module add ``spt_r11`` to the list of ``likelihoods``


Models
======

This module required ``Models`` for 

- cl_TT
- Extra-galactic foregrounds


Parameters
==========

[spt_r11].a90
-------------
[spt_r11].a150
--------------
[spt_r11].a220
--------------

    These are the calibration factors at each frequency, defined so that 
    they multiply the theory spectrum. Calibration priors
    are included in the likelihood and are 1.75%, 1.6%, and 2.4% respectively.
    
    
Plotting
========

You can use the Inspector module to plot best-fit models and residuals. 
::

    import pycosmc
    p = pycosmc.inspect('param.ini')
    p['_likelihoods']['spt_r11'].plot(p=p,delta=delta)

where ``delta`` is ``True`` if you want to plot residuals, otherwise ``False``.
"""

from spt_r11 import spt_r11
