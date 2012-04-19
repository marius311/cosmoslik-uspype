import os, sys, re
from tempfile import mkdtemp
from pycosmc.params import read_ini
from numpy import zeros, loadtxt
from pycosmc.modules import Model

class cmb_camb(Model):
    
    def init(self,p):
        pcamb = p.get('camb',{})
        if 'executable' not in pcamb: raise Exception("Expected [camb]{executable=/path/to/camb}")
        if not os.path.exists(pcamb['executable']): raise Exception("Could not find camb executable at '%s'"%pcamb['executable'])
        self.workdir = os.path.abspath(pcamb['workdir'] if 'workdir' in pcamb else mkdtemp())
        self.paramfile = os.path.join(self.workdir,'params.ini')
        self.cambini = read_ini(pcamb['defaults']) if 'defaults' in pcamb else {}
        self.cambini.update(pcamb)
        self.cambini.update(
                  {'output_root':'',
                   'total_output_file':'',
                   'lensed_total_output_file':'',
                   'lens_potential_output_file':'',
                   'scalar_output_file':os.path.join(self.workdir,'scalar'),
                   'tensor_output_file':os.path.join(self.workdir,'tensor'),
                   'lensed_output_file':os.path.join(self.workdir,'lensed'),
                   'vector_output_file':os.path.join(self.workdir,'vector'),
                   'transfer_filename(1)': os.path.join(self.workdir,'transfer'),
                   'transfer_matterpower(1)':os.path.join(self.workdir,'matterpower')
                   })
        self.outputfiles = [v for k,v in self.cambini.items() if 'output_file' in k or k in ['transfer_filename(1)','transfer_matterpower(1)']]
    
        #Read the CAMB source to find which parameters CAMB reads out of the ini file
        if pcamb.get('check_camb',False):
            self.camb_keys=set()
            for f in os.listdir(pcamb['source']):
                if f.endswith('90'):
                    with open(os.path.join(pcamb['source'],f)) as f:
                        for line in f:
                            r = re.search("Ini_Read.*File\(.*?,'(.*)'",line,re.IGNORECASE)
                            if r: self.camb_keys.add(r.group(1))
                            r = re.search("Ini_Read.*\('(.*)'",line,re.IGNORECASE)
                            if r: self.camb_keys.add(r.group(1))
            print 'Valid CAMB keys: %s'%self.camb_keys
            
    
    def get(self,p,required):
        pcamb = p.get('camb',{})
        
        Alens = p.get('Alens',1)
        self.cambini['get_scalar_cls'] = any(x in required for x in ['cl_TT','cl_TE','cl_EE','cl_BB'])
        self.cambini['get_tensor_cls'] = dotens = (p.get('r',0) != 0)
        self.cambini['get_transfer'] = dotrans = 'pk' in required
        self.cambini['do_lensing'] = dolens = (Alens != 0) 
        lmax = p['lmax']
        self.cambini['l_max_scalar'] = lmax + (100 if dolens else 0)
        lmax_tens = self.cambini['l_max_tensor'] = p.get('lmax_tensor',lmax)
        
        if not (self.cambini['get_scalar_cls'] or self.cambini['get_tensor_cls'] or self.cambini['get_transfer']): return {}
        
        #Write CAMB ini
        self.cambini.update(pcamb)
        self.cambini.update(p)
        if pcamb.get('check_camb',False): 
            print 'Setting the following CAMB parameters: %s'%{k:self.cambini[k] for k in self.cambini if isinstance(k,str) and re.sub('\([0-9]*\)','',k) in self.camb_keys}
        pcamb['ini']={}
        with open(os.path.join(self.paramfile),'w') as f:
            for k,v in self.cambini.items():
                if isinstance(v,bool): v = 'T' if v else 'F'                   
                elif isinstance(v,(float, int, str)): pass
                else: v=None
                if v!=None:
                    f.write('%s = %s\n'%(k,v))
                    pcamb['ini'][k]=v
        
        result = {}
        for x in ['TT','TE','EE','BB']: result['cl_%s'%x] = zeros(lmax)

        #Delete existing files
        for f in self.outputfiles+[os.path.join(self.workdir,'camb.out')]: 
            if os.path.exists(f): os.remove(f)

        #Call CAMB and load output files
        res = os.system('(cd %s && ./%s %s >> %s)'%(os.path.dirname(pcamb['executable']),os.path.basename(pcamb['executable']),self.paramfile,os.path.join(self.workdir,'camb.out')))
        try: pcamb['output'] = ''.join(open(os.path.join(self.workdir,'camb.out')).readlines())
        except: pcamb['output'] = None
        if res==0:
            try:
                scal = dict(zip(['l','TT','EE','TE','pp','pT'],loadtxt(self.cambini['scalar_output_file']).T))
                if dolens: lens = dict(zip(['l','TT','EE','BB','TE'],loadtxt(self.cambini['lensed_output_file']).T))
                if dotens: tens = dict(zip(['l','TT','EE','BB','TE'],loadtxt(self.cambini['tensor_output_file']).T))
                if dotrans: result['pk'] = loadtxt(self.cambini['transfer_matterpower(1)'])
            except Exception as e:
                raise Exception("Error reading CAMB output files.\n'"+str(e)+"'\nCAMB output:\n"+pcamb['output'])
        elif res!=2:
            raise Exception('CAMB returned error '+str(res)+':\n'+pcamb['output'])
        else: sys.exit()
        
        #Add scalar/lensed contribution
        for x in ['TT','TE','EE']: result['cl_%s'%x][2:lmax] += ((1-Alens)*scal[x][:lmax-2] + Alens*lens[x][:lmax-2])
        result['cl_BB'][2:lmax] += Alens*lens['BB'][:lmax-2]
        
        #Add tensor contribution
        if p.get('r',0)!=0:
            for x in ['TT','TE','EE','BB']: result['cl_%s'%x][2:lmax_tens] += tens[x][:lmax_tens-2]  
            
        return result