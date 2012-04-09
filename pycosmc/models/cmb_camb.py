import os, sys, re
from tempfile import mkdtemp
from pycosmc.ini import read_ini
from numpy import zeros, loadtxt, arange, vstack
from itertools import chain

def init(p):
    global workdir, paramfile, cambp, outputfiles, camb_keys
    
    assert 'camb' in p, "Expected camb=/path/to/camb."
    assert os.path.exists(p['camb']), "Could not find camb executable at '%s'"%p['camb']
    workdir = os.path.abspath(p['camb.workdir'] if 'camb.workdir' in p else mkdtemp())
    paramfile = os.path.join(workdir,'cambp.ini')
    cambp = read_ini(p['camb.defaults']) if 'camb.defaults' in p else {}
    cambp.update(
              {'output_root':'',
               'total_output_file':'',
               'lensed_total_output_file':'',
               'lens_potential_output_file':'',
               'scalar_output_file':os.path.join(workdir,'scalar'),
               'tensor_output_file':os.path.join(workdir,'tensor'),
               'lensed_output_file':os.path.join(workdir,'lensed'),
               'vector_output_file':os.path.join(workdir,'vector'),
               'transfer_filename(1)': os.path.join(workdir,'transfer'),
               'transfer_matterpower(1)':os.path.join(workdir,'matterpower')
               })
    outputfiles = [v for k,v in cambp.items() if 'output_file' in k or k in ['transfer_filename(1)','transfer_matterpower(1)']]

    #Read the CAMB source to find which parameters CAMB reads out of the ini file
    if p.get('check_camb',False):
        camb_keys=set()
        for f in os.listdir(p['camb.source']):
            if f.endswith('90'):
                with open(os.path.join(p['camb.source'],f)) as f:
                    for line in f:
                        r = re.search("Ini_Read.*File\(.*?,'(.*)'",line,re.IGNORECASE)
                        if r: camb_keys.add(r.group(1))
                        r = re.search("Ini_Read.*\('(.*)'",line,re.IGNORECASE)
                        if r: camb_keys.add(r.group(1))
        

def get(p,derivative=0):
        if derivative!=0: raise NotImplementedError("CAMB model can't do derivatives yet.")
        
        Alens = p.get('Alens',1)
        cambp['get_scalars'] = any(x in p['_models.get'] for x in ['cl_TT','cl_TE','cl_EE','cl_BB'])
        cambp['get_tensors'] = dotens = (p.get('r',0) != 0)
        cambp['get_transfers'] = dotrans = 'pk' in p['_models.get']
        cambp['do_lensing'] = dolens = (Alens != 0) 
        lmax = p['lmax']
        cambp['l_max_scalar'] = lmax + (100 if dolens else 0)
        lmax_tens = cambp['l_max_tensor'] = p.get('lmax_tensor',lmax)
        
        if not (cambp['get_scalars'] or cambp['get_tensors'] or cambp['get_transfers']): return {}
        
        #Write CAMB ini
        cambp.update(p)
        if p.get('check_camb'): 
            print 'Setting the following CAMB parameters: %s'%{k:p[k] for k in p if re.sub('\([0-9]*\)','',k) in camb_keys}
        p['camb.ini']={}
        with open(os.path.join(paramfile),'w') as f:
            for k,v in cambp.items():
                if isinstance(v,bool): v = 'T' if v else 'F'                   
                elif isinstance(v,(float, int, str)): pass
                else: v=None
                if v!=None:
                    f.write('%s = %s\n'%(k,v))
                    p['camb.ini'][k]=v
        
        result = {}
        for x in ['TT','TE','EE','BB']: result['cl_%s'%x] = zeros(lmax)

        #Delete existing files
        for f in outputfiles+[os.path.join(workdir,'camb.out')]: 
            if os.path.exists(f): os.remove(f)

        #Call CAMB and load output files
        res = os.system('(cd %s && %s %s >> %s)'%(os.path.dirname(p['camb']),p['camb'],paramfile,os.path.join(workdir,'camb.out')))
        try: p['camb.output'] = ''.join(open(os.path.join(workdir,'camb.out')).readlines())
        except: p['camb.output'] = None
        if res==0:
            try:
                scal = dict(zip(['l','TT','EE','TE','pp','pT'],loadtxt(cambp['scalar_output_file']).T))
                if dolens: lens = dict(zip(['l','TT','EE','BB','TE'],loadtxt(cambp['lensed_output_file']).T))
                if dotens: tens = dict(zip(['l','TT','EE','BB','TE'],loadtxt(cambp['tensor_output_file']).T))
                if dotrans: result['pk'] = loadtxt(cambp['transfer_matterpower(1)'])
            except Exception as e:
                raise Exception("Error reading CAMB output files.\n'"+str(e)+"'\nCAMB output:\n"+p['camb.output'])
        elif res!=2:
            raise Exception('CAMB returned error '+str(res)+':\n'+p['camb.output'])
        else: sys.exit()
        
        #Add scalar/lensed contribution
        for x in ['TT','TE','EE']: result['cl_%s'%x][2:lmax] += ((1-Alens)*scal[x][:lmax-2] + Alens*lens[x][:lmax-2])
        result['cl_BB'][2:lmax] += Alens*lens['BB'][:lmax-2]
        
        #Add tensor contribution
        if p.get('r',0)!=0:
            for x in ['TT','TE','EE','BB']: result['cl_%s'%x][2:lmax_tens] += tens[x][:lmax_tens-2]  
            
        return result