from cosmoslik_plugins.models.egfs import egfs
from numpy import arange, hstack, loadtxt, zeros, exp, sqrt, pi, mean, log
import os

class camspec_cibtests(egfs):
    
    eff_fr = {217: {'dust': 223.83541917989555, 'radio': 220.51977989286704, 'tsz': 222.11240613104724}, 
              143: {'dust': 144.2345587962175, 'radio': 141.5079738746441, 'tsz': 148.44314175590353},
              100: {'dust': 102.48951023720737, 'radio': 100.52206326883656, 'tsz': 99.72700223853647},
              90 : {'dust':97.9,  'radio': 95.3,  'tsz':97.6},
              150: {'dust':153.8, 'radio': 150.2, 'tsz':152.9},
              220: {'dust':219.6, 'radio': 214.1, 'tsz':218.1}}
    
    def init(self, p): 
        self.norm_ell = float(p.get(('egfs','norm_ell'),3000))
        self.norm_ell_nrun = float(p.get(('egfs','norm_ell_nrun'),3000))
        
        self.dir = os.path.dirname(os.path.abspath(__file__))
        padding = zeros(10000)
        def norm(dl,norm=self.norm_ell): return dl/dl[norm]
        def normanddl(cl,norm=self.norm_ell): return (lambda dl: dl/dl[norm])((lambda l: cl*l*(l+1)/2/pi)(arange(cl.shape[0])))
        self.tsz_template = norm(hstack([[0,0],loadtxt(os.path.join(self.dir,"camspec_templates/tsz_143_eps0.50.dat"))[:,1],padding]))
        self.ksz_template = norm(hstack([[0,0],loadtxt(os.path.join(self.dir,"camspec_templates/cl_ksz_148_trac.dat"))[:,1],padding]))
        self.cib_template = normanddl(hstack([[0,0],loadtxt(os.path.join(self.dir,"clustered_150.dat"))[:,1],padding]))

    def get_egfs(self, p, spectra, lmax, freqs, id=None, **kwargs):
        
#        import ipdb; ipdb.set_trace()
        
        
        def planck_to_spt(planck_frs,spt_frs,comp):
            """Power in Planck band to power in SPT band"""
            (pfr1,pfr2),(sfr1,sfr2) = [[self.eff_fr[f][comp] for f in frs] for frs in (planck_frs,spt_frs)]
            alpha={'dust':3.8,'radio':-0.5}[comp]
            return (sfr1*sfr2/pfr1/pfr2)**alpha/dBdT(sfr1,pfr1)/dBdT(sfr2,pfr2)

        ell = arange(lmax)/self.norm_ell

        def ellplaw(n,nrun):
            return (arange(lmax)/self.norm_ell)**(n + nrun*log(arange(lmax)/self.norm_ell))

        def cib(fr):
            return p['a_lin_%s'%fr] * self.cib_template[:lmax] + \
                   p['a_nl1_%s'%fr] * ellplaw(p['n_cib_1'],p['nrun_cib_1']) + \
                   p['a_nl2_%s'%fr] * ellplaw(p['n_cib_2'],p['nrun_cib_2'])

        def tsz(fr1,fr2):
            return p['a_tsz'] * self.tsz_template[:lmax] * tszdep(fr1,fr2, p['tsz_norm_fr'])
        
        
        p = p['egfs']
        
        comps={}
        comps['ksz'] = p['a_ksz'] * self.ksz_template[:lmax] 
        
        if id=='planck':
        
            freqs = tuple(fr if isinstance(fr,dict) else {'tsz':fr} for fr in freqs)
            
            frlbl = tuple([n for n,(l,u) in {100:(80,120),
                                             143:(130,160), 
                                             217:(200,240)}.items() 
                           if l<mean(fr.values())<u][0] for fr in freqs)
            
            comps['tsz'] = tsz(freqs[0]['tsz'],freqs[1]['tsz'])
            
            if frlbl==(100,100):
                comps['ps_100'] = p['a_ps_100']*ell**2
                comps['cib_100'] = cib('100')
            elif frlbl==(143,143):
                comps['ps_143'] = p['a_ps_143']*ell**2
                comps['cib_143'] = cib('143')
                comps['tsz_cib'] = - 2 * p['xi'] * sqrt(tsz(143,143)*cib(143))
            elif frlbl==(217,217):
                comps['ps_217'] = p['a_ps_217']*ell**2
                comps['cib_217'] = cib('217')
            elif tuple(sorted(frlbl))==(143,217):
                comps['ps_143_217'] = p['r_ps']*sqrt(p['a_ps_143']*p['a_ps_217'])*ell**2
                comps['cib_143_217'] = p['r_cib']*sqrt(cib('143')*cib('217'))
                comps['tsz_cib'] = - p['xi'] * sqrt(tsz(143,143)*cib(217))                
            
        elif id.startswith('spt'):
            
            fr1, fr2 = freqs
            frlbl = tuple([n for n,(l,u) in {90:(80,120),
                                 150:(130,160), 
                                 220:(200,240)}.items() 
               if l<mean(fr.values())<u][0] for fr in freqs)

            
            comps['dgpo'] = p['a_dgpo_spt'] * ell**2 * (fr1['dust']*fr2['dust']/p['dgpo_spt_norm_fr']**2)**p['alpha_dgpo_spt'] / dBdT(fr1['dust'],p['dgpo_spt_norm_fr']) / dBdT(fr2['dust'],p['dgpo_spt_norm_fr'])
            comps['radio'] = p['a_radio_spt'] * ell**2 * (fr1['radio']*fr2['radio']/p['radio_spt_norm_fr']**2)**p['alpha_radio_spt'] / dBdT(fr1['radio'],p['radio_spt_norm_fr']) / dBdT(fr2['radio'],p['radio_spt_norm_fr'])
            
            comps['tsz'] = tsz(fr1['tsz'],fr2['tsz'])

            if frlbl==(150,150):
                comps['dgcl'] = cib('143') * planck_to_spt((143,143), frlbl, 'dust')
            elif frlbl==(150,220):
                comps['dgcl'] = p['r_cib'] * sqrt(cib('143')*cib('217')) * planck_to_spt((143,217), frlbl, 'dust')
            elif frlbl==(220,220):
                comps['dgcl'] = cib('217') * planck_to_spt((217,217), frlbl, 'dust')
            elif frlbl==(90,90):
                comps['dgcl'] = cib('100') * planck_to_spt((100,100), frlbl, 'dust')
            elif frlbl==(90,220):
                comps['dgcl'] = p['r_cib_90_220'] * sqrt(cib('100')*cib('217')) * planck_to_spt((100,217), frlbl, 'dust')
            elif frlbl==(90,150):
                comps['dgcl'] = p['r_cib_90_150'] * sqrt(cib('100')*cib('143')) * planck_to_spt((100,143), frlbl, 'dust')
                
        return comps
        
def tszdep(fr1,fr2,fr0):
    """The tSZ frequency dependence."""
    t1,t2,t0 = map(lambda fr: (lambda x0: x0*(exp(x0)+1)/(exp(x0)-1) - 4)(fr/56.78),[fr1,fr2,fr0])
    return t1*t2/t0**2
        
def dBdT(fr1,fr0):
    """ dB/dT at T_CMB """
    dBdT,dBdT0 = map((lambda fr: (lambda x0: x0**4 * exp(x0) / (exp(x0)-1)**2)(fr/57.78)),[fr1,fr0])
    return dBdT/dBdT0          
