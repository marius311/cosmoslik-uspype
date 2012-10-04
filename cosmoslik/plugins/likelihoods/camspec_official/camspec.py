from cosmoslik.plugins import Likelihood, SubprocessExtension
from numpy import arange, pi, hstack

class camspec_official(Likelihood):
    """
    
    ======================
    The CAMspec Likelihood
    ====================== 

    Authors: George Efstathiou, Steven Gratton
    Plugin by: Marius Millea

    A thin wrapper of the CAMspec Fortran likelihood.
    
    This plugin is not compatible with CosmoSlik extra-galactic foregrounds. 
    It uses it's own built-in foreground model controlled by a set of parameters
    in the [camspec] section of the ini. 
    
    
    Parameters:
    -----------
    
    
    [camspec]{
    
        #initialization parameters
        like_file = 
        sz143_file = 
        tszxcib_file = 
        ksz_file =
        beam_file = 
        
        #foreground parameters
        a_ps_100 =
        a_ps_143 =
        a_ps_217 = 
        a_cib_143 = 
        a_cib_217 = 
        a_sz = 
        a_ksz = 
        r_ps =
        r_cib =
        xi =
        cal0 =
        cal1 =
        cal2 =

        #beam eigen-mode amplitudes
        
        [spec_1]{
            mode_1 = 1
            mode_2 = 1
            mode_3 = 1
            ...
        }
        ...
        
    }
    
    
    """
    
    def init(self, p):
        self.pycamspec = SubprocessExtension('pycamspec',globals())
        self.pycamspec.like_init(**{x:p['camspec',x] 
                                    for x in ['like_file', 
                                              'sz143_file', 
                                              'tszxcib_file', 
                                              'ksz_file', 
                                              'beam_file']})
    
    def lnl(self, p, model):
        
        def tocl(dl): return hstack([[0],dl[1:]/arange(1,len(dl))/(arange(1,len(dl))+1)])
        
        return self.pycamspec.calc_like(cell_cmb=tocl(model['cl_TT']),
                                        beam_coeffs=[[p.get(('camspec','spec_%i'%i,'mode_%i'%j),0)
                                                      for j in range(1,self.pycamspec.num_modes_per_beam+1)]
                                                     for i in range(1,self.pycamspec.beam_nspec+1)],
                                        **{x.lower():p['camspec',x]
                                           for x in ['A_ps_100','A_ps_143','A_ps_217',
                                                     'A_cib_143','A_cib_217','A_sz','A_ksz',
                                                     'r_ps','r_cib','xi',
                                                     'cal0', 'cal1', 'cal2']})
        
