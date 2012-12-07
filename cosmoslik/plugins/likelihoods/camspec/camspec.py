from cosmoslik.plugins import Likelihood, SubprocessExtension
from numpy import arange, hstack

class camspec(Likelihood):
    """
    
    ======================
    The CAMspec Likelihood
    ====================== 
    
    Authors: George Efstathiou, Steven Gratton
    Plugin by: Marius Millea
    
    This plugin is fully compatible with CosmoSlik extra-galactic foregrounds.    
    
    
    Parameters:
    -----------
    
    
    [camspec]{
    
        #initialization parameters
        like_file = 
        sz143_file = 
        tszxcib_file = 
        ksz_file =
        beam_file = 
        
        #effective frequencies and flux cuts to be used for foreground modeling
        eff_fr = {100: {'dust':143, 'radio':143, 'tsz':143}, ...}
        fluxcut = {100: 200, ...}
        
        
        #calibration parameters
        cal0 =
        cal1 =
        cal2 =

        #beam eigen-mode amplitudes
        [spec_0]{
            mode_0 = 1
            mode_1 = 1
            mode_2 = 1
            ...
        }
        ...
        
    }
    
    
    """
    
    def init(self, p):
        from pycamspec import pycamspec
        self.pycamspec = pycamspec
#        self.pycamspec = SubprocessExtension('pycamspec',globals())
        self.pycamspec.like_init(**{x:p['camspec',x] 
                                    for x in ['like_file', 
                                              'sz143_file', 
                                              'tszxcib_file', 
                                              'ksz_file', 
                                              'beam_file']})
        self.lmax = max(self.pycamspec.lmaxx)
        self.freqs = p['camspec','eff_fr']
        self.fluxcut = p['camspec','fluxcut']
        
    
    def get_required_models(self,p):
        return ['cl_TT','egfs']
    
    
    def get_cl_model(self, p, model=None):
        if model is None: model = p['_model']
        clmodel = {'cmb':model['cl_TT']}
        clmodel.update({(fr1,fr2):model['egfs']('cl_TT',
                                                lmax=self.lmax, 
                                                freqs=(self.freqs[fr1],self.freqs[fr2]),
                                                fluxcut=min(self.fluxcut[fr1],self.fluxcut[fr2]))
                        for fr1,fr2 in [('100','100'),('143','143'),('217','217'),('143','217')]})
        return clmodel
    
    
    def lnl(self, p, model):
        def tocl(dl): return hstack([[0],dl[1:]/arange(1,len(dl))/(arange(1,len(dl))+1)])[:self.lmax]
        cl_model = self.get_cl_model(p, model)
        return self.pycamspec.calc_like(cell_cmb=tocl(cl_model['cmb']),
                                        cell_fg100=tocl(cl_model[('100','100')]),
                                        cell_fg143=tocl(cl_model[('143','143')]),
                                        cell_fg217=tocl(cl_model[('217','217')]),
                                        cell_fg143x217=tocl(cl_model[('143','217')]),
                                        beam_coeffs=[[p.get(('camspec','spec_%i'%i,'mode_%i'%j),0)
                                                      for j in range(self.pycamspec.num_modes_per_beam)]
                                                     for i in range(self.pycamspec.beam_nspec)],
                                        **{x:p.get(('camspec',x),1) for x in ['cal0', 'cal1', 'cal2']})
        