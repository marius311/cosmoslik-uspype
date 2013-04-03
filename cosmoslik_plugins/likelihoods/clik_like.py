from cosmoslik.plugins import Likelihood
from numpy import hstack, zeros

class clik_like(Likelihood):
    """
    
    [clik]{
    
        files = {'highL':'path...'}
    
        [highl]{
            nuisance...
        }
    
    }
    
    TODO: implement nuisance parameters
    
    """
    
    def init(self, p):
        import clik
        self.files = p.get(('clik','files'),{})
        self.cliks = {k:clik.clik(v) for k,v in self.files.items()}
        
    def get_required_models(self,p):
        return {'cl_%s'%x  for clik in self.cliks.values() for x, lmax in zip(['TT','EE','BB','TE','TB','EB'],clik.get_lmax()) if lmax!=-1}
            
    def lnl(self, p, model):
        lnl = 0
        for k,clik in self.cliks.items():
            lnl += clik(hstack(model.get('cl_%s'%x,zeros(lmax+1))[:lmax+1] 
                               for x, lmax in zip(['TT','EE','BB','TE','TB','EB'],clik.get_lmax()) 
                               if lmax!=-1))[0]
        
        return lnl