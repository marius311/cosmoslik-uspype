import mspec as M
from numpy import dot, arange, diag
from scipy.linalg import cho_solve, cho_factor
from pycosmc.modules import Likelihood


class mspec_lnl(Likelihood):
    
    
    def get_required_models(self,p):
        return ['cl_TT']


    def init(self,p):
        if 'mspec' not in p: raise Exception('Expected an [mspec] section in the ini file.')
        
        self.mp = M.read_Mspec_ini(p['mspec'])
        
        self.signal = M.load_signal(self.mp).dl()
        
        #All the things we the per-frequency C_ell signal
        self.process_signal = lambda s: s.lincombo(self.mp['cleaning']) \
                                         .rescaled(self.mp.get('rescale',1)) \
                                         .sliced(self.mp['binning'](slice(*self.mp["lrange"])))
        
        self.processed_signal = self.process_signal(self.signal)
        (self.signal_matrix_spec, self.signal_matrix_cov) = self.processed_signal.get_as_matrix(ell_blocks=True)
        
        self.signal_matrix_cov = cho_factor(self.signal_matrix_cov)
        self.signal_matrix_spec = self.signal_matrix_spec[:,1]
        
        self.fluxcut = self.mp['fluxcut']
        self.eff_fr = self.mp['eff_fr']
        self.lmax = self.mp["lrange"][1]
        
    
    def lnl(self,p,model):
        
        #Build an Mspec PowerSpectra object for all the required frequencies 
        #which holds CMB + foreground C_ell's
        model_sig = M.PowerSpectra(ells=arange(self.lmax))
        for fr1,fr2 in self.signal.get_spectra():
            cl = model['cl_TT'][:self.lmax]
            cl += model['egfs']('cl_TT',
                               fluxcut=min(self.fluxcut[fr1],self.fluxcut[fr2]),
                               freqs=(self.eff_fr[fr1],self.eff_fr[fr2]),
                               lmax=self.lmax)
            model_sig[(fr1,fr2)] = model_sig[(fr2,fr1)] = cl
        
        
        #Diagnostic plotting
        if p.get('diagnostic',False):
            from matplotlib.pyplot import ion, draw, cla, yscale, ylim
            ion()
            cla()
            self.processed_signal.plot()
            self.process_signal(model_sig.binned(self.mp['binning'])).plot()
            yscale('log')
            ylim(10,6e3)
            draw()

        #Do all the same binning, slicing, etc... that we do to the signal
        model_sig = self.process_signal(model_sig.binned(self.mp['binning'])).get_as_matrix(ell_blocks=True).spec[:,1]
        
        #Compute the likelihood  
        dcl = model_sig - self.signal_matrix_spec
        return dot(dcl,cho_solve(self.signal_matrix_cov,dcl))/2
