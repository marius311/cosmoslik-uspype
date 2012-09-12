from cosmoslik.plugins import Likelihood
from numpy import loadtxt, sum, array, pi
import os

class planck_cib(Likelihood):
    
    def init(self,p):
        
        self.dir = os.path.dirname(os.path.abspath(__file__))
        
        self.dat545 = loadtxt(os.path.join(self.dir,'dat','CIB_avg_545_545_pk_K_WithProjBeam.dat'))
        self.ells = self.dat545[:,0]
        self.dat545[:,[3,4]] = (self.dat545[:,[3,4]].T * self.ells**2/2/pi * 1e12).T
        self.binslices = [slice(*x) for x in array(self.dat545[:,[1,2]],dtype=int)]
        
    def get_required_models(self,p):
        return ['egfs']
    
    
    def lnl(self,p,model):
        
        return sum((self.get_cl_model(p)- self.dat545[:,3])**2/2/self.dat545[:,4]**2)
    
    def get_cl_model(self, p):
        clcib = p['_model']['egfs']('cl_TT',
                                    freqs = ({k:545 for k in ['dust','radio','tsz']},)*2,
                                    fluxcut = 200,
                                    lmax = 3500)
        return [clcib[x].mean() for x in self.binslices]


    def plot(self,p,ax=None,residuals=True):
        if ax is None:
            from matplotlib.pyplot import figure
            ax = figure().add_subplot(111)

        if residuals:
            ax.errorbar(self.ells,self.dat545[:,3] - self.get_cl_model(p),yerr=self.dat545[:,4],ls='',marker='.')
        else:
            ax.plot(self.ells, self.get_cl_model(p))
            ax.errorbar(self.ells,self.dat545[:,3],yerr=self.dat545[:,4],ls='',marker='.')