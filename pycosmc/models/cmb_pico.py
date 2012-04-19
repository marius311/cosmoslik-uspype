import pypico
import pycosmc.models.cmb_camb
from pycosmc.modules import Model

class cmb_pico(Model):
    
    num_pico = 0
    num_camb = 0
    
    def init(self,p):
        try: datafile = p['pico','datafile']
        except KeyError: raise Exception("Please specify [pico]{datafile = ... }")
        else: self.pico = pypico.loadpico(*datafile)
            
        self.camb = pycosmc.models.cmb_camb.cmb_camb()
        self.camb.init(p)
        
    def get(self,p,required):
        if (self.num_pico+self.num_camb)%10==0 and p.get('pico_verbose',False): print 'PICO=%i CAMB=%i'%(self.num_pico,self.num_camb)
        try: 
            r = self.pico.get(outputs=required,**p)
            self.num_pico+=1
            return r
        except pypico.CantUsePICO as e: 
            self.num_camb+=1
            print e
            print 'Calling CAMB...'
            return self.camb.get(p,required)
