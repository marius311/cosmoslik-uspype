from ast import literal_eval
from numpy import array
import re

def read_ini(path):
   
    params = {}
    sampled_params = []
   
    with open(path) as f: lines=[l for l in [re.sub("#.*","",l).strip() for l in f.readlines()] if len(l)>0]
    for line in lines:
        tokens = [t.strip() for t in line.split("=")]
        if (len(tokens)!=2):
            raise SyntaxError("Error parsing "+file+". Expected one \"=\" in line '"+line+"'")
        elif (tokens[0][0] in ["$","*"]):
            raise SyntaxError("Error parsing "+file+". Key can't start with "+tokens[0][0]+". '"+line+"'")
        elif tokens[1]!='':
            k, v = tokens[0], tokens[1]
            if (type(v)==str):
                r = re.search("({0})\s\[\s?({0})\s({0})\s({0})\s?\]".format("[0-9.eE+-]+?"),v)
                if (r!=None):
                    v=float(r.groups()[0])
                    params["*"+k]=map(float,r.groups())
                    sampled_params.append(k)
                else:
                    try: 
                        v = literal_eval(v)
                        if isinstance(v, list): v=array(v)
                    except: pass
                
            params[k] = v
        
    params["_sampled"] = sampled_params 
    params["$OUTPUT"] = sampled_params + params.get("derived","").split()
    
    return params

def add_sampled_param(p, name, value, min, max, width, output=True):
    p[name]=value
    p['*'+name]=[value,min,max,width]
    p['_sampled'].append(name)
    if output: p['$OUTPUT'].append(name)
    
    