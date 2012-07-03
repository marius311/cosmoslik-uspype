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
        
        self.processed_signal = self.process_signal(self.signal)
        (self.signal_matrix_spec, self.signal_matrix_cov) = self.processed_signal.get_as_matrix(ell_blocks=True)
        
        self.signal_matrix_cov = cho_factor(self.signal_matrix_cov)
        self.signal_matrix_spec = self.signal_matrix_spec[:,1]
        
        self.fluxcut = self.mp['fluxcut']
        self.eff_fr = self.mp['eff_fr']
        self.lmax = self.mp["lrange"][1]
        
    def process_signal(self,s):
        """All the thing we do to the signal after loading it in."""
        if 'cleaning' in self.mp: s=s.lincombo(self.mp['cleaning'])
        return s.rescaled(self.mp.get('rescale',1)) \
                .sliced(self.mp['binning'](slice(*self.mp["lrange"])))

    def get_cl_model(self,p,model):
        """ 
        Build an Mspec PowerSpectra object which holds CMB + foreground C_ell's
        for all the required frequencies 
        """
        model_sig = M.PowerSpectra(ells=arange(self.lmax))
        for fr1,fr2 in self.signal.get_spectra():
            cl = model['cl_TT'][:self.lmax].copy()
            cl += model['egfs']('cl_TT',
                               fluxcut=min(self.fluxcut[fr1],self.fluxcut[fr2]),
                               freqs=(self.eff_fr[fr1],self.eff_fr[fr2]),
                               lmax=self.lmax)
            model_sig[(fr1,fr2)] = model_sig[(fr2,fr1)] = cl

        return self.process_signal(model_sig.binned(self.mp['binning']))


    def plot(self,ax,key,cl=None,p=None):
        if cl==None: cl=self.get_cl_model(p, p['_model'])
        if key=='cl_TT':
            self.processed_signal.plot(ax=ax)
            cl.plot(ax=ax)
        elif key=='delta_cl_TT':
            self.processed_signal.diffed(cl.spectra.values()[0]).plot(ax=ax)
            ax.plot([cl.ells[0],cl.ells[-1]],[0]*2)


    
    def lnl(self,p,model):
        
        cl_model = self.get_cl_model(p, model)
        
        cl_model_matrix = cl_model.get_as_matrix(ell_blocks=True).spec[:,1]
        
        #Diagnostic plotting
        if p.get('diagnostic',False):
            from matplotlib.pyplot import ion, draw, cla, yscale, ylim
            ion()
            cla()
#            ax.set_yscale('log')
#            ax.set_ylim(10,6e3)
            draw()
        
        #Compute the likelihood  
        dcl = cl_model_matrix - self.signal_matrix_spec
        return dot(dcl,cho_solve(self.signal_matrix_cov,dcl))/2
