from xs_lookup import LIB,find_cell,xs_at,load_nuclide,sigma_at
import os
import h5py
import numpy as np

def show_mt301(nuc,T):
    path=os.path.join(LIB,nuc+".h5")
    with h5py.File(path,"r")as f:
        root=f[nuc]
        names=list(root["reactions"].keys())
        for name in sorted(names):
            print("name=%s"%root["reactions"][name])
        print("reaction_301"in names)
        print("reaction_901"in names)
        if "reaction_301"in names:
            redundant=float(root["reactions"]["reaction_301"].attrs["redundant"])
            print("redundant=%d"%redundant)
            thi=float(root["reactions"]["reaction_301"][T]["xs"].attrs["threshold_idx"])
            print("threshold_idx=%d"%thi)
            xs=root["reactions"]["reaction_301"][T]["xs"][:]
            print("len(xs)=%d"%len(xs))
            energy=root["energy"][T][:]
            print("len(energy network)=%d"%len(energy))

def load_heating(nuc,T):
    path=os.path.join(LIB,nuc+".h5")
    with h5py.File(path,"r")as f:
        root=f[nuc]
        E=root["energy"][T][:]
        xs_301=root["reactions"]["reaction_301"][T]["xs"][:]
        thi_301=int(root["reactions"]["reaction_301"][T]["xs"].attrs["threshold_idx"])
                # [Claude assert] 抓：读错了温度组 / 阈值下标没对上 → heat 和能量网格错位，查出来的值全是隔壁格点的
        assert len(xs_301) + thi_301 == len(E), \
            "%s 的 MT=301 有 %d 个点、阈值下标 %d，加起来应等于能量网格的 %d 个点" % (nuc, len(xs_301), thi_301, len(E))
        heat=np.zeros(len(E))
        heat[thi_301:]=xs_301
    return E,heat    
            
if __name__=="__main__":
    nucs=["H1","O16","B10"]
    T="294K"
    # show_mt301(nuc,T)
    E_TEST=[0.0253,1.0e6]
    for nuc in sorted(nucs):
        print("nuc is %s"%nuc)
        E,heat=load_heating(nuc,T)
        E_net,mts,XS,awr,Q=load_nuclide(nuc,T)
        for energy in E_TEST:
            lo=find_cell(E,energy)
            xs_energy=xs_at(E,heat,energy,lo)
            xs,_=sigma_at(E_net,XS,energy)
            for i in range(len(xs)):
                print("MT=%d,xs=%.6f"%(mts[i],xs[i]))
            print("energy=%e,xs_energy(MT=301)=%e"%(energy,xs_energy))
        print("len(E)=%d,len(heat)=%d"%(len(E),len(heat)))
        print("heat[first]=%.6f,heat[last]=%.6f"%(heat[0],heat[-1]))
