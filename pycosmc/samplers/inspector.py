from pycosmc.modules import Sampler
from pycosmc.post import load_chain, Chains

class inspector(Sampler):
    
    def init(self, p):
        c = load_chain(p.pop('output_file'))
        if isinstance(c,Chains): c=c.join()
        self.bestfit = c.best_fit()
        
    def sample(self, x, lnl, p):
        x = [self.bestfit['.'.join(k)] for k in p.get_all_sampled()]
        lnl, p = lnl(x, p)
        yield x, 1, lnl, p
