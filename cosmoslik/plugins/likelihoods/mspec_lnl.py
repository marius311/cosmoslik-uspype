import mspec as M
from numpy import dot, arange, diag
from scipy.linalg import cho_solve, cho_factor, cholesky
from cosmoslik.plugins import Likelihood
from itertools import combinations_with_replacement

class mspec_lnl(Likelihood):
    
    
    def get_required_models(self,p):
        return ['cl_TT']

    def init(self,p):
        if 'mspec' not in p: raise Exception('Expected an [mspec] section in the ini file.')
        
        self.mp = M.read_Mspec_ini(p['mspec'])
        
        self.signal = M.load_signal(self.mp)#.dl()
        
        self.processed_signal = self.process_signal(self.signal)

        self.lrange = self.mp['lrange']
        if isinstance(self.lrange,tuple): 
            self.lrange = {k:self.lrange for k in self.processed_signal.get_spectra()}
        
        (self.signal_matrix_spec, self.signal_matrix_cov) = self.processed_signal.get_as_matrix(lrange=self.lrange)
        
        self.signal_matrix_cov = cholesky(self.signal_matrix_cov), False
        self.signal_matrix_spec = self.signal_matrix_spec[:,1]
        
        self.fluxcut = self.mp['fluxcut']
        self.eff_fr = self.mp['eff_fr']
        self.lmax = max([u for (_,u) in self.lrange.values()])
        

        
    def process_signal(self,s):
        """All the thing we do to the signal after loading it in."""
        if 'cleaning' in self.mp: s=s.lincombo(self.mp['cleaning'])
        s = s.rescaled(self.mp.get('rescale',1))
        return s
    
    
    def get_cl_model(self,p,model=None):
        """ 
        Build an Mspec PowerSpectra object which holds CMB + foreground C_ell's
        for all the required frequencies 
        """
        if model is None: model = p['_model']
        model_sig = M.PowerSpectra(ells=arange(self.lmax))
        for fr1,fr2 in self.signal.get_spectra():
            cl = model['cl_TT'][:self.lmax].copy()
            cl += model['egfs']('cl_TT',
                               fluxcut=min(self.fluxcut[fr1],self.fluxcut[fr2]),
                               freqs=(self.eff_fr[fr1],self.eff_fr[fr2]),
                               lmax=self.lmax)
            model_sig[(fr1,fr2)] = model_sig[(fr2,fr1)] = cl

        return self.process_signal(model_sig.binned(self.mp['binning']))


    def plot(self,
             fig=None,
             cl=None, 
             p=None, 
             show_comps=False,
             yscale='log',
             ylim=None,
             residuals=False):
        
        if cl==None: cl=self.get_cl_model(p, p['_model'])
        if fig==None: 
            from matplotlib.pyplot import figure
            fig=figure()
            
        n=len(self.processed_signal.get_maps())

        fig.set_size_inches(6*n,6*n/1.6)
        fig.subplots_adjust(hspace=0,wspace=0)
            
        for ((i,fri),(j,frj)) in combinations_with_replacement(enumerate(self.processed_signal.get_maps()),2):
            ax=fig.add_subplot(n,n,n*j+i+1)
            if residuals:
                self.processed_signal.diffed(cl[(fri,frj)]).plot(ax=ax,which=[(fri,frj)],c='k')
                ax.plot([cl.ells[0],cl.ells[-1]],[0]*2)
            else:
                self.processed_signal.sliced(self.processed_signal.binning(slice(*self.lrange[(fri,frj)]))).plot(ax=ax,which=[(fri,frj)],c='k')
                cl.plot(ax=ax,which=[(fri,frj)],c='k')
                if show_comps:
                    ax.plot(p['_model']['cl_TT'],c='b')
                    p['_model']['egfs']('cl_TT',
                                        fluxcut=min(self.fluxcut[fri],self.fluxcut[frj]),
                                        freqs=(self.eff_fr[fri],self.eff_fr[frj]),
                                        lmax=self.lmax,
                                        plot=True,
                                        ax=ax)
                    
                ax.set_ylim(*(ylim or ((0,6999) if yscale=='linear' else (11,9999))))
                ax.set_yscale(yscale)
                ax.set_xlim(2,self.lmax-1)
                
                if n==1:
                    ax.set_title('%sx%s'%(fri,fri))
                else:
                    if i==0: ax.set_ylabel(frj,size=16)
                    else: ax.set_yticklabels([])
                    if j==n-1: ax.set_xlabel(fri,size=16)
                    else: ax.set_xticklabels([])



    
    def lnl(self,p,model):
        
        cl_model = self.get_cl_model(p, model)
        
        cl_model_matrix = cl_model.get_as_matrix(lrange=self.lrange).spec[:,1]
        
        #Compute the likelihood  
        dcl = cl_model_matrix - self.signal_matrix_spec
        return dot(dcl,cho_solve(self.signal_matrix_cov,dcl))/2
