from pycosmc.modules import Sampler
from pycosmc.post import load_chain, Chains

class inspector(Sampler):
    """
    
    [inspector]
        axes = ['ax1','ax2']
        [ax1]
            content = ('spt_k11','cl_TT')
            ylim
            xlim
        [ax2]
            content = ('spt_k11','cl_TT')
            ylim
            xlim

    
    """
    
    def init(self, p):
        c = load_chain(p.pop('output_file'))
        if isinstance(c,Chains): c=c.join()
        self.bestfit = c.best_fit()
        
    def sample(self, x, lnl, p):
        x = [self.bestfit['.'.join(k)] for k in p.get_all_sampled()]
        lnl, p = lnl(x, p)
        yield x, 1, lnl, p
