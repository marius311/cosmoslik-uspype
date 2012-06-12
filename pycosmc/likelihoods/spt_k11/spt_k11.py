from numpy import array, fromstring, loadtxt, dot, arange, diag
from scipy.linalg import cho_factor, cho_solve
from pycosmc.modules import Likelihood
import os

class spt_k11(Likelihood):

    def lnl(self, p, model):
        #Get CMB + foreground model
        cl = model['cl_TT'][:self.lmax] + \
             model['egfs']('cl_TT', lmax=self.lmax, freqs=(self.freq,self.freq), fluxcut=self.fluxcut)
        
        #Apply window functions
        cl = array([dot(cl[self.windowrange],w) for w in self.windows])
            
        if p.get('diagnostic',False):
            from matplotlib.pyplot import ion, errorbar, plot, draw, cla, yscale, ylim
            ion()
            cla()
            errorbar(self.ells,self.spec,yerr=diag(self.sigma[0]),fmt='.',label='SPT K11')
            plot(self.ells,cl)
            yscale('log')
            ylim(10,6e3)
            draw()

            
        #Apply windows and calculate likelihood
        dcl = self.spec-cl
        return dot(dcl,cho_solve(self.sigma, dcl))/2

    
    def get_required_models(self, p):
        return ['cl_TT', 'egfs']

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
        self.lmax = self.windowrange.stop
        self.ells = array([dot(arange(10000)[self.windowrange],w) for w in self.windows])

        self.freq = {'dust':154, 'radio': 151, 'tsz':153}
        self.fluxcut = 50
        
    




