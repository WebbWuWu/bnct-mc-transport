# -*- coding: utf-8 -*-
# ============================================================
#  photon.py  --  M4：γ 光子输运（光子截面、Klein–Nishina、光子栈）
#  文件头（import、常量）Claude 2026-11-08 写；函数你写。
#  框架：项目文档\M4框架_γ光子输运_2026-11-08.md
#  第 ① 段用法：项目文档\load_photon用法_2026-11-08.md
#  光子数据：lib80x_hdf5\photon\H.h5、O.h5、B.h5（ENDF/B-VIII.0 光原子库，11/7 explore\make_photon_lib.py 转的；和 OpenMC 同一份）
# ============================================================
import os, math
import h5py
import numpy as np

_LIB_DIR = r"E:\BNCT_GPU\lib80x_hdf5" if os.name == "nt" else "/mnt/e/BNCT_GPU/lib80x_hdf5"   # Windows / WSL 都能 import
PHOTON_DIR = os.path.join(_LIB_DIR, "photon")

PHOTON_ELEMS = ["H", "O", "B"]                       # 水 + 硼里用到的元素（光子只认元素，不认同位素）
NUC2ELEM = {"H1": "H", "O16": "O", "B10": "B"}       # nuclides 里的 "name" → 光子库的元素名

def loglog0(x,y,e):
    lo=np.searchsorted(x,e,side="right")-1
    if y[lo]<=0.0 or y[lo+1]<=0.0:
        return 0.0
    else:
        frac=(math.log(e)-math.log(x[lo]))/(math.log(x[lo+1])-math.log(x[lo]))
        lny=math.log(y[lo])+frac*(math.log(y[lo+1])-math.log(y[lo]))
        return math.exp(lny)
    
def load_photon():
    photon={}
    for elem in PHOTON_ELEMS:
        path=os.path.join(PHOTON_DIR,elem+".h5")
        with h5py.File(path,"r")as f:
            g=f[elem]
            photon[elem]={}
            photon[elem]["E"]=g["energy"][()]
            photon[elem]["coh"]=g["coherent/xs"][()]
            photon[elem]["incoh"]=g["incoherent/xs"][()]
            photon[elem]["pe"]=g["photoelectric/xs"][()]
            photon[elem]["pair"]=g["pair_production_nuclear/xs"][()]+g["pair_production_electron/xs"][()]
    return photon

def photon_sigma(photon,nuclides,E):
    sig={"coh":0.0,"incoh":0.0,"pe":0.0,"pair":0.0}
    for nuc in nuclides:
        elem=NUC2ELEM[nuc["name"]]
        dens=nuc["dens"]
        for k in list(sig.keys()):
            sig[k]+=loglog0(photon[elem]["E"],photon[elem][k],E)*(1e-24)*dens
    sig["tot"]=sig["coh"]+sig["incoh"]+sig["pe"]+sig["pair"]
    return sig
        

if __name__=="__main__":
    from lesson11_water import kur_table,KEEP_FACTOR,build_water
    pho=load_photon()
    for elm in PHOTON_ELEMS:
        print(elm,len(pho[elm]["E"]),pho[elm]["E"][0],pho[elm]["E"][-1],len(pho[elm]["pair"]))
    edges,shares=kur_table()
    nuclides,order=build_water(KEEP_FACTOR*edges[-1])
    for E in [2.2233e6,1e6,0.4778e6]:
        sig=photon_sigma(pho,nuclides,E)
        print(E,sig["tot"],sig["pair"]/sig["tot"])
