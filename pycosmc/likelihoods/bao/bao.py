from scipy.integrate import quad 
from numpy import *
from numpy.linalg import inv
from pycosmc.modules import Likelihood


class bao(Likelihood):
    """
    Beutler et al. (2011) give rs/Dv = 0.336 \pm 0.015  at z=0.106
    Padmanabhan et al. (2012) give Dv/rs = 8.88 \pm 0.17 at z=0.35
    BOSS (Anderson et al. 2012) report Dv/rs = 13.67 \pm 0.22 at z = 0.57
    Wigglez report rs/Dv = 0.0726 \pm 0.0034 at z = 0.6
    """

    def init(self,p):
        self.zvec=array([0.106,0.35,0.57,0.6])
        self.Dv_over_rs_data=array([2.98,8.88,13.67,13.77])
        self.Dv_over_rs_cov=diag([0.13,0.17,0.22,0.65])**2 #diagonal, but could change with new data      

    def lnl(self,p,model):
        dx = self.Dv_over_rs_data - [Dv_over_rs(z=z,**p) for z in self.zvec]
        return dot(dx,dot(inv(self.Dv_over_rs_cov),dx))/2
        

def Dv_over_rs(ombh2, omch2, H0, z_drag, z, **kwargs):
    rs = r_s(z=z_drag,**arguments(besides=['z']))
    return (D_A(**arguments())**2*z/Hubble(**arguments()))**(1/3.)/rs





"""
TODO: These belong in some kind of utility module
"""

def arguments(ascend=0, besides=None):
    """Returns a dictionary of all named arguments passed to this function."""
    from inspect import getargvalues, stack
    kwname, args = getargvalues(stack()[1+ascend][0])[-2:]
    args.update(args.pop(kwname, []))
    if besides!=None: 
        for b in besides:
            if b in args: args.pop(b)
    return args

    
def kw_vectorize(func):
    """Like numpy.vectorize but allows keyword arguments."""
    import inspect
    vfunc = vectorize(func)
    def myfunc(**kwargs):
        return vfunc(*[kwargs[k] for k in inspect.getargspec(func).args],**kwargs)
    return myfunc




"""
TODO: All the things below here belong in a Model 
"""

taufactor = 2.30952e-5  #combination of constants (Thomson cross section, G, baryon mass) important for optical depth calc (units=Mpc^{-1})
rhoxOverOmegaxh2 = 8.09906e-11 #eV^4
G = 1.63995e2  #G in eV^(-4)Mpc^(-2)
EightPiGOver3= 8*pi/3.*G  #in eV^(-4)Mpc^(-2)
#For neutrinos and photons
KelvinToeV=8.6173324e-5
Tgamma0=2.7255*KelvinToeV  
rhogamma0=pi**2/15.*Tgamma0**4
OneHundredKilometersPerSecPerMpcOverSpeedofLightTimesMpc = 3.33565e-4


def r_s(z, ombh2, **kwargs):
    Rovera=3.*ombh2*rhoxOverOmegaxh2/(4.*rhogamma0)
    allkw = arguments(besides=['z'])
    return quad(lambda zp: 1/Hubble(zp,**allkw) / sqrt(3*(1+Rovera/(1+zp))),z,inf)[0]
    
def D_A(omkh2,**kwargs):
    """
    dist is comoving proper distance (calculated by distance prop)
    returns comoving angular-diameter distance in Mpc
    """
    K=-omkh2*(OneHundredKilometersPerSecPerMpcOverSpeedofLightTimesMpc)**2
    dist = D_prop(**arguments())
    if K<0: return 1./sqrt(-K)*sin(dist*sqrt(-K))
    elif K > 0: return 1./sqrt(K)*sinh(dist*sqrt(K))
    elif K==0: return dist


def D_prop(z,ommh2,omkh2,omvh2,mnu,Nnu_massive,Nnu_massless,**kw):
    """returns proper comoving distance to redshift z in Mpc"""
    allkw = arguments(besides=['z'])
    return quad(lambda zp: 1/Hubble(z=zp,**allkw),0,z)[0]


def Hubble(z,ommh2,omkh2,omvh2,mnu,Nnu_massive,Nnu_massless,**kw):
    """
    Returns H in Mpc^{-1}
    m is mass of single neutrino species in eV
    Nmass is number of massive species (assumed to each have same mass)
    Neff is number of massless species
    """
    ae=1.e-5
    Te=Tgamma0/ae*(4./11.)**(1./3.) #scale photon temp, then convert to neutrino temp
    Nphotoneff=1.+Nnu_massless*7./8.*(4./11.)**(4./3.) #effective number of photon species for massless radiation
    
    a = 1./(1.+z)
    rhonu = 0 if mnu==0 else rhoredshift(a,ae,Te,mnu/(1.*Nnu_massive))  #contribution from one massive neutrino species

    return sqrt(EightPiGOver3*( rhoxOverOmegaxh2*( ommh2*a**(-3) + omkh2*a**(-2) + omvh2 ) + Nphotoneff*rhogamma0*a**(-4)+ Nnu_massive*rhonu ))
    
    
def rhoredshift(a,ae,Te,m):
    """
    ae is early scale factor when still relativistic,
    Te is temperature at that time in eV
    m  is mass in eV
    a  is scale factor at epoch for which rho is desired
       output:  rho (in units of ev^4)
    """
#    if ae/a*Te < 0.01*m: rho=1./(8.*pi**(3./2.))*exp(-m/Te)*m*(2*m*Te)**(1.5)*(ae/a)**3
    integral = quad(lambda x: x**2*sqrt(x**2*(ae/a)**2+m**2)/(exp(sqrt(x**2+m**2)/Te)+1.0),0.0, inf)[0]
    return 1./(2.*pi**2)*(ae/a)**3*integral


def taub(ombh2, omch2, Yp, ommh2, omkh2, omvh2, mnu, Nnu_massive, Nnu_massless, Xe, **kw):
    """ The integrand for the taub integration. """
    allkw = arguments(besides=['z'])
    
    Rovera=3.*ombh2*rhoxOverOmegaxh2/(4*rhogamma0)  #where R = 3\rho_b/(4\rho_\gamma)
    
    @vectorize
    def taub(z,zstart=0):
        return taufactor*(1.-Yp) * quad(lambda zp: Xe(zp)*(1.+zp)**2/(1.+Rovera/(1+zp))/Hubble(z=zp,**allkw),zstart,z,full_output=1)[0]
    
    return taub

                   
