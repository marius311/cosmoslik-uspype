from numpy import array, loadtxt, dot, arange, diag, hstack, zeros
from scipy.linalg import cho_factor, cho_solve
from pycosmc.modules import Likelihood
import os

class spt_r11(Likelihood):

    def lnl(self, p, model):
            
        cl = self.get_cl_model(model)
            
        if p.get('diagnostic',False):
            from matplotlib.pyplot import ion, figure, draw, cla
            ion()
            cla()
            ax = figure(0).add_subplot(111)
            self.plot(ax,'cl_TT',cl)
            ax.set_yscale('log')
            ax.set_ylim(10,6e3)
            draw()
            
        cl_vector = hstack([cl[spec_name] for spec_name in self.spec_names])

        
        dcl = self.spec_vector - cl_vector
        return dot(dcl,cho_solve(self.cov, dcl))/2


    def plot(self,ax,key,cl=None,p=None):
        if cl==None: cl = self.get_cl_model(p['_model'])
        if key=='cl_TT':
            for spec_name, c in zip(self.spec_names,'bgrckm'):
                ax.errorbar(self.ells[spec_name],self.spec[spec_name],yerr=self.sigmas[spec_name],c=c,fmt='.',label='SPT K11')
                ax.plot(self.ells[spec_name],cl[spec_name],c=c)


    def get_cl_model(self,model):
        def get_cl(fr1,fr2): return hstack([model['cl_TT'],zeros(10000)])[:self.lmax] + model['egfs']('cl_TT', lmax=self.lmax, freqs=(self.freq[fr1],self.freq[fr2]), fluxcut=self.fluxcut)
        def apply_windows(cl,windows): return array([dot(cl[self.windowrange],w[:,1]) for w in windows])
        return {spec_name:apply_windows(get_cl(*spec_name),windows) for (spec_name, windows) in self.windows.items()}
    
    def to_matrix(self,spec):
        return hstack([spec[spec_name] for spec_name in self.spec_names])
        
    def get_required_models(self, p):
        return ['cl_TT', 'egfs']
    


    def init(self, p):
        
        self.datadir     = os.path.join(os.path.dirname(__file__),'spt_multif_0809')
        
        self.spec_names  = [('90','90'),('90','150'),('90','220'),('150','150'),('150','220'),('220','220')]
        self.spec_vector = loadtxt(os.path.join(self.datadir,"spectrum_90_90x150_90x220_150_150x220_220.txt"))[:,1]
        self.spec        = dict(zip(self.spec_names,self.spec_vector.reshape(6,15)))
        self.cov         = cho_factor(loadtxt(os.path.join(self.datadir,"covariance_90_90x150_90x220_150_150x220_220.txt")).reshape(90,90))
        self.sigmas      = dict(zip(self.spec_names,diag(self.cov[0]).reshape(6,15)))
        self.windows     = dict(zip(self.spec_names,(lambda w: w.reshape(6,15,w.shape[1],2))(array([loadtxt(os.path.join(self.datadir,'window_%i'%i)) for i in range(1,91)]))))
        self.windowrange = (lambda w: slice(w[0,0],w[-1,0]+1))(self.windows.values()[0][0])
        self.ells        = {frs: array([dot(arange(self.windowrange.start,self.windowrange.stop),w[:,1]) for w in windows]) for frs,windows in self.windows.items()} 
        self.lmax        = self.windowrange.stop
        self.freq        = {'90' : {'dust':97.9,  'radio': 95.3,  'tsz':97.6},
                            '150': {'dust':153.8, 'radio': 150.2, 'tsz':152.9},
                            '220': {'dust':219.6, 'radio': 214.1, 'tsz':218.1}}
        self.fluxcut     = 6.4

        


