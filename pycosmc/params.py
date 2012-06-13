from ast import literal_eval
from numpy import array, inf
import re
from collections import OrderedDict

class SectionDict(dict):
    """ Dictionary class for PyCosMC with a few extra convenience methods. """
    
    def section(self,key):
        """ Returns the dictionary for a subsection, or {} if that subsection is not present. """
        return self.get(key,{})
    
    def update(self,other):
        for k,v in (other.iteritems() if isinstance(other,dict) else other): self[k]=v
    
    def __getitem__(self, key):
        if type(key)==tuple and len(key)==1: key=key[0]
        if type(key)==tuple:
            try: return super(SectionDict,self).__getitem__(key[0])[key[1:]]
            except KeyError: raise KeyError(key)
        else: 
            return super(SectionDict,self).__getitem__(key)
        
    def __setitem__(self, key, value):
        if type(key)==tuple and len(key)==1: key=key[0]
        if type(key)==tuple:
            if key[0] not in self: self[key[0]] = SectionDict()
            super(SectionDict,self).__getitem__(key[0])[key[1:]]=value
        else: 
            return super(SectionDict,self).__setitem__(key,value)
                
    def get(self, key, default=None):
        try: return self[key]
        except KeyError: return default
        
    def setdefault(self, key, default=None):
        if default==None: default = SectionDict()
        try: return self[key]
        except KeyError:
            self[key]=default
            return default

    def __contains__(self, key):
        try: self[key]; return True
        except KeyError: return False
    
    def add_sampled_param(self, name, value, min=-inf, max=inf, width=1):
        """ Add a sampled parameter. """
        self[name]=value
        self.setdefault('_sampled',{})[name] = [value,min,max,width]
        
    def get_all_sampled(self):
        """ Recursively get sampled params in all sub-sections. """
        sampled = [((k,),v) for k,v in self.get('_sampled',{}).iteritems()]
        for k, v in self.iteritems():
            if isinstance(v,SectionDict): 
                sampled += [((k,)+k2, v2) for k2, v2 in v.get_all_sampled().iteritems()] 
        return OrderedDict(sorted(sampled))
    
    def copy(self):
        return SectionDict(self)


def read_ini(f):   
    """ Load parameters from ini """ 
    
    if isinstance(f,str): f=open(f)
    f = enumerate(f,1)
    
    def _read_ini(f,subgroup=False):
        d=SectionDict()
        for i,line in f:
            if not line: break
            line = re.sub('#.*','',line).strip()
            if line=='}':
                if subgroup: return d
                else: raise ValueError("Unexpected '}' at line %i"%i)
            elif line!='':
                r = re.match('\[(.*)\]{$',line)
                if r!=None: 
                    d[r.group(1)]=_read_ini(f,True)
                    continue

                r = re.match('include (.*)',line)
                if r!=None:
                    d.update(read_ini(r.group(1)))
                    continue

                r = re.match('\s*(.*?)\s*=\s*(.*?)\s*$',line)
                if r!=None:
                    d[r.group(1)]=r.group(2)
                    continue
                    
                raise ValueError("Could not parse line %i: '%s'"%(i,line))
                                
        if subgroup: raise ValueError("Expected a '}'")
        else: return d
    
    return _read_ini(f, False)


def eval_values(p):
    """ Try to evaluate values as Python expressions. """
    for k,v in p.items():
        if isinstance(v,dict): eval_values(v)
        else:
            try: 
                v = literal_eval(v)
                if isinstance(v, list): v=array(v)
                p[k]=v
            except: pass
        
        
    
def process_parameters(params):
    """ Process parameters. """
   
    params['_sampled'] = {}
   
    for k,v in params.iteritems():
        if (type(v)==str):
            r = re.search("({0})\s\[\s?({0})\s({0})\s({0})\s?\]".format("[0-9.eE+-]+?"),v)
            if (r!=None):
                params.add_sampled_param(k,*map(float,r.groups()))
            else:
                try: 
                    v = literal_eval(v)
                    if isinstance(v, list): v=array(v)
                    params[k] = v
                except: pass
        elif isinstance(v,SectionDict):
            process_parameters(v)
            
