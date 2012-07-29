from cosmoslik.modules import Likelihood
from numpy import loadtxt, sum, array, arange
import os

class planck_lens(Likelihood):
    
    
    
    
    def init(self,p):
        
        self.dir = os.path.dirname(os.path.abspath(__file__))
        
        self.lrange = array(p.get(('planck_lens','lrange'),[100,501]))
        
        #these files start at ell=2
        self.clpp_bestfit = loadtxt(os.path.join(self.dir,'bestfit_scalCls.dat'))[slice(*(self.lrange-2)),4]
        self.plm_norm = loadtxt(os.path.join(self.dir,'plm_norm.dat'))[slice(*(self.lrange-2)),1]
        self.plm_nhl = loadtxt(os.path.join(self.dir,'dat_plm_nhl.dat'))[slice(*(self.lrange-2)),1]
        
        self.denominator = sum(self.clpp_bestfit**2 * (2.*arange(*self.lrange)+1.) / 2. / (self.plm_norm**2 * self.plm_nhl)**2)
        
        
    def get_required_models(self,p):
        return ['cl_pp']
    
    
    def lnl(self,p,model):
        """
        clpp is an array containing containing l^4 C_l^{\phi \phi} (standard CAMB format) at every l starting at l = 0 and up
        to some undetermined lmax, with lmax > 500.  It will only be accessed for 99 < l < 501.
        Adat is the measured value of what Duncan calls A_{100}^{500}.  Amod is the model value for this number.  
        We caclulate Amod by averaging clpp/clpp_bestfit.  We then compare to Adat to get lnL
        """
        A_data = .992  # Can get Planck value in Jan 23, 2012 email from DH 
        sigma_A_data = 0.061 #(or slide 16 of  http://wiki.planck.fr/index.php/Meetings/2012-01-18?action=download&upname=jplensing.pdf)
    
        return ((A_data - self.get_A_model(model['cl_pp']))/sigma_A_data)**2/2.

    
    def get_A_model(self, clpp):
        return sum(self.clpp_bestfit * clpp[slice(*self.lrange)] * (2.*arange(*self.lrange)+1.) / 2. / (self.plm_norm**2 * self.plm_nhl)**2) / self.denominator
    
