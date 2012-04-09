import pypico
import pycosmc.models.cmb_camb as C

num_pico = 0
num_camb = 0

def init(p):
    global pico
    if ("pico_datafile" in p):
        print "Initializing PICO..."
        pico = pypico.loadpico(*p["pico_datafile"])
    else:
        raise Exception("Please specify pico_datafile.")
    C.init(p)
    
def get(p,derivative=0):
    if derivative!=0: raise NotImplementedError('CMB model derivative not implemented yet.')
    global num_pico, num_camb
    if (num_pico+num_camb)%10==0 and p.get('pico_verbose',False): print 'PICO=%i CAMB=%i'%(num_pico,num_camb)
    try: 
        r = pico.get(**p)
        num_pico+=1
        return r
    except pypico.CantUsePICO as e: 
        num_camb+=1
        print e
        print 'Calling CAMB...'
        return C.get(p,derivative=derivative)
