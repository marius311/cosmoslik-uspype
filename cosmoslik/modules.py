class Likelihood:
    """ A cosmoslik likelihood module. """
    
    def init(self,p):
        """ Do any initialization for the likelihood here. """
        pass

    def lnl(self,p,model):
        """ Return negative log-likelihood. """
        raise NotImplementedError('Implement Likelihood.lnl.')

    def get_extra_params(self,p):
        """ Return a list of nuisance parameter names needed by this likelihood. """
        return {}
        
    def get_required_models(self,p):
        """ Return a list of model components needed by this likelihood, e.g. 'cl_TT' """
        return []
        
        
        
class Model:
    """ A cosmoslik model module. """
    
    def init(self,p):
        """ Do any initialization for the model here. """
        pass

    def get(self,p,required):
        """ Get the model. """
        raise NotImplementedError('Implement Model.get.')
    
    

class Sampler:
    """ A cosmoslik sampler. """
    
    def init(self,p):
        """ Do any initialization for the sampler here. """
        pass
        
    def sample(self, x, lnl, p):
        """ Return a generator which yields samples. """
        raise NotImplementedError('Implement Sampler.sample')
    
    
class Deriver:
    
    def init(self,p):
        """ Do any initialization for the deriver here. """
        
    def add_derived(self,p):
        """ Add derived parameters. """
        raise NotImplementedError('Implement Deriver.add_derived')
