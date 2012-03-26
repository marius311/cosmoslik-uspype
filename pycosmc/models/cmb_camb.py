import os
from tempfile import mkdtemp
from subprocess import check_output
from ini import read_ini
from numpy import zeros, loadtxt, array

def init(p):
    assert os.path.exists(p['camb']), "Could not find camb executable at '%s'"%p['camb']
    if 'camb.workdir' not in p: p['camb.workdir'] = mkdtemp()
    p['camb.paramfile'] = os.path.join(p['camb.workdir'],'params.ini')
    p['camb.params'] = read_ini(p['camb.defaults']) if 'camb.defaults' in p else {}
    p['camb.params'].update(
              {'output_root':'',
               'total_output_file':'',
               'lensed_total_output_file':'',
               'lens_potential_output_file':'',
               'scalar_output_file':os.path.join(p['camb.workdir'],'scalar'),
               'tensor_output_file':os.path.join(p['camb.workdir'],'tensor'),
               'lensed_output_file':os.path.join(p['camb.workdir'],'lensed'),
               'vector_output_file':os.path.join(p['camb.workdir'],'vector'),
               'transfer_filename(1)': os.path.join(p['camb.workdir'],'transfer'),
               'transfer_matterpower(1)':os.path.join(p['camb.workdir'],'matterpower')
               })


def get(p,derivative=0):
        result = {}
        
        Alens = p.get('A_lens',1)
        p['camb.params']['do_lensing'] = dolens = (Alens != 0) 
        p['camb.params']['get_tensors'] = dotens = (p.get('r',0) != 0)
        p['camb.params']['l_max_scalar'] = p['lmax'] + 100 if dolens else 0
        p['camb.params']['l_max_tensor'] = p['lmax_tensor'] if 'lmax_tensor' in p else p['lmax']
        
        #Write CAMB ini
        with open(os.path.join(p['camb.paramfile']),'w') as f:
            for k,v in dict(p['camb.params'],**p).items():
                if isinstance(v,bool): f.write('%s = %s\n'%(k,'T' if v else 'F'))
                elif isinstance(v,(float, int, str)): f.write('%s = %s\n'%(k,v))

        #Call CAMB
        check_output([p['camb'],p['camb.paramfile']],cwd=os.path.dirname(p['camb']))
    
        #Load CAMB output files
        scal = dict(zip(['l','TT','EE','TE','pp','pT'],loadtxt(p['camb.params']['scalar_output_file']).T))
        if dolens: lens = dict(zip(['l','TT','EE','BB','TE'],loadtxt(p['camb.params']['lensed_output_file']).T))
        if dotens: tens = dict(zip(['l','TT','EE','BB','TE'],loadtxt(p['camb.params']['tensor_output_file']).T))
        nls, nlt = min(p['lmax'],len((lens if dolens else scal)['l'])), min(p['lmax'],len(tens['l']) if dotens else 0)

        for x in ['TT','TE','EE','BB']: result['cl_%s'%x] = zeros(p['lmax']+2)
        
        #Add scalar/lensed contribution
        for x in ['TT','TE','EE']: result['cl_%s'%x][2:nls+2] += (1-Alens)*scal[x][:nls] + Alens*lens[x][:nls]
        result['cl_BB'][2:nls+2] += Alens*lens['BB'][:nls]
        
        #Add tensor contribution
        if p.get('r',0)!=0:
            for x in ['TT','TE','EE','BB']: result['cl_%s'%x][2:nlt+2] += tens[x][:nlt]  
            
        return result
        

#p={'lmax':1000,'camb.workdir':os.path.abspath(''),'camb':'/home/marius/workspace/camb/camb/camb','camb.defaults':'/home/marius/workspace/camb/camb/params.ini'}
#init(p)
#get(p)
#from matplotlib.pyplot import * 
#for k,v in p['model'].items(): 
#    semilogy(v,label=k)
#legend(loc='lower right')
#show()