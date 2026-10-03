import os
import h5py
import numpy as np

LIB=r"E:\BNCT_GPU\lib80x_hdf5"
NUC="B10"
T="294K"
E_TEST=0.0253
BARN=1.0e-24

def load_nuclide(nuc,T):
    path=os.path.join(LIB,nuc+".h5")
    with h5py.File(path,"r")as f:
        root=f[nuc]
        E=root["energy"][T][:]
        awr=float(root.attrs["atomic_weight_ratio"])
        names=list(root["reactions"].keys())
        
        mts=[]
        XS=[]
        
        for name in sorted(names):
            redundant=root["reactions"][name].attrs["redundant"]
            if redundant==1:
                continue
            xs_ds=root["reactions"][name][T]["xs"]
            xs=xs_ds[:]
            thr=int(xs_ds.attrs["threshold_idx"])
            s=np.zeros(len(E))
            s[thr:]=xs
            XS.append(s)
            mts.append(int(root["reactions"][name].attrs["mt"]))
    return E,mts,XS,awr

def xs_at(E,s,e,lo):
    i=lo+1
    frac=(e-E[lo])/(E[i]-E[lo])
    return s[lo]+frac*(s[i]-s[lo])

def find_cell(E,e):
    if e==E[-1]:
        lo=len(E)-2
    else:
        i=np.searchsorted(E,e,side="right")
        lo=i-1
    return lo

def sigma_at(E,XS,e):
    lo=find_cell(E,e)
    sig=[]
    for j in range(len(XS)):
        sig.append(xs_at(E,XS[j],e,lo))
    total=sum(sig)
    return sig,total

if __name__ == "__main__":
    E, mts, XS, awr = load_nuclide(NUC, T)
    i2   = mts.index(2)
    i800 = mts.index(800)
    i801 = mts.index(801)
    i51  = mts.index(51)
    lo = find_cell(E, E_TEST)
    print("elastic  :", xs_at(E, XS[i2],   E_TEST,lo))
    print("n,a0     :", xs_at(E, XS[i800], E_TEST,lo))
    print("n,a1     :", xs_at(E, XS[i801], E_TEST,lo))
    print("n,a total:", xs_at(E, XS[i800], E_TEST,lo) + xs_at(E, XS[i801], E_TEST,lo))
    print("inelastic:", xs_at(E, XS[i51],  E_TEST,lo))
