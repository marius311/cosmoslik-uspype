from cosmoslik.plugins import Likelihood
from numpy import hstack, zeros, arange, pi, inf

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
        tot_lnl = 0
        
        for k,clik in self.cliks.items():
            lnl=-clik(hstack(tocl(model.get('cl_%s'%x,zeros(lmax+1))[:lmax+1])
                             for x, lmax in zip(['TT','EE','BB','TE','TB','EB'],clik.get_lmax()) 
                             if lmax!=-1))
            if lnl==0: return inf
            else: tot_lnl += lnl
            
        return lnl
    
def tocl(dl): 
    return hstack([zeros(2),dl[2:]/arange(2,dl.size)/(arange(2,dl.size)+1)*2*pi])
    