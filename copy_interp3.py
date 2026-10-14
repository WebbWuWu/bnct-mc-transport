import math
import numpy as np

NA   = 6.02214076e23
M_N  = 1.00866491588
RHO  = 1.06

def find(E,e):
    n=len(E)
    
    if e>=E[n-1]:
        return n-2
    if e<=E[0]:
        return 0
    lo=0
    hi=n-1
    step=0
    while True:
        step+=1
        mid=(hi-lo)//2+lo
        
        if e>E[mid]:
            lo=mid
        if e<E[mid]:
            hi=mid
        if e==E[mid]:
            return mid
        if hi-lo==1:
            return lo

def interp_loglog(E_tab,S_tab,e):
    k=int(np.searchsorted(E_tab,e,side="right"))-1
    k=min(max(k,0),len(E_tab)-2)
    p=np.log(S_tab[k+1]/S_tab[k])/np.log(E_tab[k+1]/E_tab[k])
    return S_tab[k]*(e/E_tab[k])**p

def number_density(w,awr):
    M=awr*M_N
    return w*RHO*NA/M

if __name__ == "__main__":
    E = [0.37, 1.9, 4.6, 12.5, 31.0, 88.0]
    S = [41.0, 18.2, 11.7, 7.1, 4.9, 3.3]
    for e in [0.21, 0.37, 7.3, 12.5, 88.0, 140.0]:
        print("find  e=%-6s -> k=%d" % (e, find(E, e)))
    for e in [4.6, 7.3, 12.5]:
        print("loglog e=%-5s -> %.10f" % (e, interp_loglog(E, S, e)))
    print("N = %.6e cm^-3" % number_density(0.708, 15.8575))
