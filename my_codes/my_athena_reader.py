
import numpy as np

def read_tab(fname):
    
    with open(fname, "r") as fn:
        
        data = {}
        
        h1,h2 = fn.readlines()[:2]
        
        for word in h1.split():
            if "time" in word:
                time = float(word.split("=")[-1])
                data['time'] = time
            if "cycle" in word:
                cycle = int(word.split("=")[-1])
                data['cycle'] = cycle
                
        vrs = h2.split()[1:]
        
    fldata = np.genfromtxt(fname,skip_header=2)
    
    fldata = fldata.T
    
    for i,vr in enumerate(vrs):
        data[vr] = fldata[i]
    
    return data
