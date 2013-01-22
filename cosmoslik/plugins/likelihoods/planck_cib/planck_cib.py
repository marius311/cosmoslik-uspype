from cosmoslik.plugins import Likelihood
from numpy import loadtxt, sum, array, pi
import os

class planck_cib(Likelihood):
    
    def init(self,p):
        
        self.dir = os.path.dirname(os.path.abspath(__file__))
        self.frs = p.get(('planck_cib','freqs'),('545','545'))
        self.dat = loadtxt(os.path.join(self.dir,'dat','CIB_avg_%s_pk_K_WithProjBeam.dat'%('_'.join(self.frs))))
        self.ells = self.dat[:,0]
        self.dat[:,[3,4]] = (self.dat[:,[3,4]].T * self.ells**2/2/pi * 1e12).T
        self.binslices = [slice(*x) for x in array(self.dat[:,[1,2]],dtype=int)]
        
    def get_required_models(self,p):
        return ['egfs']
    
    
    def lnl(self,p,model):
        
        return sum((self.get_cl_model(p)- self.dat[:,3])**2/2/self.dat[:,4]**2)
    
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
            ax.errorbar(self.ells,self.dat[:,3] - self.get_cl_model(p),yerr=self.dat[:,4],ls='',marker='.')
        else:
            ax.plot(self.ells, self.get_cl_model(p))
            ax.errorbar(self.ells,self.dat[:,3],yerr=self.dat[:,4],ls='',marker='.')
