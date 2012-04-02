import os
from tempfile import mkdtemp
from subprocess import check_output
from pycosmc.ini import read_ini
from numpy import zeros, loadtxt

def init(p):
    global workdir, paramfile, params, outputfiles
    
    assert 'camb' in p, "Expected camb=/path/to/camb."
    assert os.path.exists(p['camb']), "Could not find camb executable at '%s'"%p['camb']
    workdir = os.path.abspath(p['camb.workdir'] if 'camb.workdir' in p else mkdtemp())
    paramfile = os.path.join(workdir,'params.ini')
    params = read_ini(p['camb.defaults']) if 'camb.defaults' in p else {}
    params.update(
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
    outputfiles = [v for k,v in params.items() if 'output_file' in k or k in ['transfer_filename(1)','transfer_matterpower(1)']]

def get(p,derivative=0):
        if derivative!=0: raise NotImplementedError("CAMB model can't do derivatives yet.")
        
        Alens = p.get('A_lens',1)
        params['get_scalars'] = any(x in p['_models.get'] for x in ['cl_TT','cl_TE','cl_EE','cl_BB'])
        params['get_tensors'] = dotens = (p.get('r',0) != 0)
        params['get_transfers'] = dotrans = 'pk' in p['_models.get']
        params['do_lensing'] = dolens = (Alens != 0) 
        params['l_max_scalar'] = p['lmax'] + 100 if dolens else 0
        params['l_max_tensor'] = p['lmax_tensor'] if 'lmax_tensor' in p else p['lmax']
        
        if not (params['get_scalars'] or params['get_tensors'] or params['get_transfers']): return {}
        
        #Write CAMB ini
        with open(os.path.join(paramfile),'w') as f:
            for k,v in dict(params,**p).items():
                if isinstance(v,bool): f.write('%s = %s\n'%(k,'T' if v else 'F'))
                elif isinstance(v,(float, int, str)): f.write('%s = %s\n'%(k,v))

        #Delete existing files
        for f in outputfiles: 
            if os.path.exists(f): os.remove(f)

    
        result = {}
        for x in ['TT','TE','EE','BB']: result['cl_%s'%x] = zeros(p['lmax']+2)

        #Call CAMB and load output files
        try:
            output = check_output([p['camb'],paramfile],cwd=os.path.dirname(p['camb']))
            scal = dict(zip(['l','TT','EE','TE','pp','pT'],loadtxt(params['scalar_output_file']).T))
            if dolens: lens = dict(zip(['l','TT','EE','BB','TE'],loadtxt(params['lensed_output_file']).T))
            if dotens: tens = dict(zip(['l','TT','EE','BB','TE'],loadtxt(params['tensor_output_file']).T))
            if dotrans: result['pk'] = loadtxt(params['transfer_matterpower(1)'])
        except Exception:
            raise Exception('CAMB error\n'+output)
            
        nls, nlt = min(p['lmax'],len((lens if dolens else scal)['l'])), min(p['lmax'],len(tens['l']) if dotens else 0)
        
        #Add scalar/lensed contribution
        for x in ['TT','TE','EE']: result['cl_%s'%x][2:nls+2] += (1-Alens)*scal[x][:nls] + Alens*lens[x][:nls]
        result['cl_BB'][2:nls+2] += Alens*lens['BB'][:nls]
        
        #Add tensor contribution
        if p.get('r',0)!=0:
            for x in ['TT','TE','EE','BB']: result['cl_%s'%x][2:nlt+2] += tens[x][:nlt]  
            
        #Get matter power-spectrum 
        
        return result