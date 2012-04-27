import os, sys, re
from pycosmc.params import read_ini
from numpy import zeros, loadtxt
from pycosmc.modules import Model
from cStringIO import StringIO
from subprocess import Popen, PIPE

class camb(Model):
    
    def init(self,p):
        pcamb = p.get('camb',{})
        self.cambdir = os.path.abspath(os.path.join(os.path.dirname(__file__),'camb'))
        pcamb.setdefault('executable',os.path.join(self.cambdir,'camb'))
        if not os.path.exists(pcamb['executable']): raise Exception("Could not find camb executable at '%s'"%pcamb['executable'])
        self.cambdefs = read_ini(pcamb['defaults']) if 'defaults' in pcamb else {}
        self.cambdefs.update(pcamb)
        self.cambdefs.update(
                  {'output_root':'',
                   'total_output_file':'',
                   'lensed_total_output_file':'',
                   'lens_potential_output_file':'',
                   'scalar_output_file':'scal',
                   'tensor_output_file':'tens',
                   'lensed_output_file':'lens',
                   'vector_output_file':'vec',
                   'transfer_filename(1)': 'trans',
                   'transfer_matterpower(1)':'mpk'
                   })
        self.cambproc = Popen(["./camb"],cwd=self.cambdir,stdin=PIPE,stdout=PIPE,stderr=PIPE)
    
        #Read the CAMB source to find which parameters CAMB reads out of the ini file
        if pcamb.get('check_camb',False):
            self.camb_keys=set()
            for f in os.listdir(pcamb['source']):
                if f.endswith('90'):
                    with open(os.path.join(pcamb.get('source',self.cambdir),f)) as f:
                        for line in f:
                            r = re.search("Ini_Read.*File\(.*?,'(.*)'",line,re.IGNORECASE)
                            if r: self.camb_keys.add(r.group(1))
                            r = re.search("Ini_Read.*\('(.*)'",line,re.IGNORECASE)
                            if r: self.camb_keys.add(r.group(1))
            print 'Valid CAMB keys: %s'%self.camb_keys
            
    
    def get(self,p,required):
        pcamb = p.get('camb',{})
        
        Alens = p.get('Alens',1)
        cambini = pcamb['ini'] = dict(self.cambdefs)
        cambini['get_scalar_cls'] = doscal = any(x in required for x in ['cl_TT','cl_TE','cl_EE','cl_BB'])
        cambini['get_tensor_cls'] = dotens = (p.get('r',0) != 0)
        cambini['get_transfer'] = dotrans = 'pk' in required
        cambini['do_lensing'] = dolens = (doscal and Alens != 0)
        docl = doscal or dolens or dotens 
        lmax = pcamb['lmax']
        cambini['l_max_scalar'] = lmax + (100 if dolens else 0)
        lmax_tens = cambini['l_max_tensor'] = p.get('lmax_tensor',lmax)
        
#        if not (doscal or dolens or dotens or dotrans): return {}
        
        #Write CAMB ini
        cambini.update(p)
        if pcamb.get('check_camb',False): 
            print 'Setting the following CAMB parameters: %s'%{k:self.cambini[k] for k in self.cambini if isinstance(k,str) and re.sub('\([0-9]*\)','',k) in self.camb_keys}
        
        for k,v in cambini.items():
            if isinstance(v,bool): pcamb['ini'][k] = 'T' if v else 'F'                   
            elif isinstance(v,(float, int, str)): pcamb['ini'][k] = v
            else: cambini.pop(k)
        
        result = {}
        cambini = '\n'.join('%s = %s'%(k,v) for k,v in pcamb['ini'].items()) + '\n$\n'

        #Send params to CAMB and read output
        self.cambproc.stdin.write(cambini)
        output = read_camb_output(self.cambproc.stdout, self.cambproc.stderr)
        
        if doscal: scal = dict(zip(['l','TT','EE','TE','pp','pT'],output['scal'].T))
        if dolens: lens = dict(zip(['l','TT','EE','BB','TE'],output['lens'].T))
        if dotens: tens = dict(zip(['l','TT','EE','BB','TE'],output['tens'].T))
        if dotrans: result['mpk'] = output['mpk']
        
        #Combine cl contributions
        if docl:
            for x in ['TT','TE','EE','BB']: 
                if 'cl_%s'%x in required:
                    result['cl_%s'%x] = zeros(lmax)
                    if doscal or dolens: 
                        result['cl_%s'%x][2:lmax] += (((1-Alens)*scal[x][:lmax-2] if x!='BB' and doscal else 0)) + (Alens*lens[x][:lmax-2] if dolens else 0)
                    if dotens:
                        result['cl_%s'%x][2:lmax_tens] += tens[x][:lmax_tens-2]
        
        #TODO: figure out where to put z_drag
        p['z_drag'] = output['z_drag']
        
        return result        


def read_camb_output(stdout, stderr):
    out, outdict = '', {}
    cur = None
    while True:
        line = stdout.readline()
        out += (line+'\n')
        if line=='': raise Exception('CAMB error.\n%s%s'%(out,'\n'.join(stderr.readlines())))
        if line.strip()=='$': break
        r = re.match('\[(.*)\]',line.strip())
        if r!=None: 
            cur=r.group(1)
        else:
            matches = list(re.finditer('\s*(.+?)\s*=\s*(.+?)(\s|$)',line))
            if len(matches)>0:
                for m in matches: outdict[m.group(1)]=m.group(2)
            else:
                outdict[cur] = outdict.get(cur,'') + line + '\n'
                    
    for k, v in outdict.iteritems():
        try: outdict[k]=loadtxt(StringIO(v))
        except: pass
    
    return outdict

