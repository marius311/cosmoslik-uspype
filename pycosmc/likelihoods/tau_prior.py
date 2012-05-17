from pycosmc.modules import Likelihood

class tau_prior(Likelihood):
    
    def lnl(self, p, model):
        return (p['tau'] - 0.0851)/(2*0.014**2)
    
    
