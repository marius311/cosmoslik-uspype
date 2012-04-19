from numpy import array, fromstring, loadtxt, dot, arange, diag
from scipy.linalg import cho_factor, cho_solve
from pycosmc.modules import Likelihood
import os

class spt_k11(Likelihood):

    def get_required_models(self, p):
        return ['cl_TT','pk']


    def get_extra_params(self,p):
        return {'Aps':[100, 0, 300, 10]}
    
    
    def init(self, p):
        
        self.datadir = os.path.join(os.path.dirname(__file__),'bandpowers')
        
        #Load spectrum and covariance
        with open(os.path.join(self.datadir,'Spectrum_spt20082009.newdat')) as f:
            while 'TT' not in f.readline(): pass
            self.spec=array([fromstring(f.readline(),sep=' ')[1] for _ in range(47)])
            self.sigma=cho_factor(array([fromstring(f.readline(),sep=' ') for _ in range(94)])[47:])
            
        #Load windows
        self.windows = [loadtxt(os.path.join(self.datadir,'windows','window_0809','window_%i'%i))[:,1] for i in range(1,48)]
        self.windowrange = (lambda x: slice(min(x),max(x)+1))(loadtxt(os.path.join(self.datadir,'windows','window_0809','window_1'))[:,0])
        self.ells = array([dot(arange(10000)[self.windowrange],w) for w in self.windows])
        
    
    
    def lnl(self, p, model):
        #Get CMB + foreground model
        cl = model['cl_TT'].copy()
        if len(cl)<self.windowrange.stop: raise Exception("SPT K11 likelihood needs C_ell's to ell=%i"%self.windowrange.stop)
        cl += p['spt_k11','Aps']*(arange(len(cl))/3000.)**2 #Hacked PS term until egfs module
        #if 'fgs' in model: cl += model['fgs'](eff_fr=eff_fr,fluxcut=fluxcut)
        cl = array([dot(cl[self.windowrange],w) for w in self.windows])
        
        #Apply windows and calculate likelihood
        dcl = self.spec-cl
        
        if p.get('diagnostic',False):
            from matplotlib.pyplot import ion, errorbar, plot, draw, cla, yscale, ylim
            ion()
            cla()
            errorbar(self.ells,self.spec,yerr=diag(self.sigma[0]),fmt='.',label='SPT K11')
            plot(self.ells,cl)
            yscale('log')
            ylim(10,6e3)
            draw()
            
        return dot(dcl,cho_solve(self.sigma, dcl))/2
