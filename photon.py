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
from rng import isotropic_direction
from lib_path import LIB_DIR
from thermal import rotate_dir
_LIB_DIR = LIB_DIR                                   # 11/10：库路径统一写在 lib_path.py（值没变）
PHOTON_DIR = os.path.join(_LIB_DIR, "photon")

PHOTON_ELEMS = ["H", "O", "B"]                       # 水 + 硼里用到的元素（光子只认元素，不认同位素）
NUC2ELEM = {"H1": "H", "O16": "O", "B10": "B"}       # nuclides 里的 "name" → 光子库的元素名
E_G_H=2.2233e6
E_G_B=0.4778e6
HC=12398.419843320026
ME_C2=0.51099895e6
E_CUT_G=1.0e3

def record_site(p,xs,hist):
    sig_cap=xs["gH_coe"]+xs["gB_coe"]
    y=sig_cap/xs["Sigma_t"]
    share=xs["gH_coe"]/sig_cap
    dic={
        "z":p["z"],
        "wt":p["wt"],
        "y":y,
        "share":share,
    }
    hist["gsites"].append(dic)

def track_photon(g, hist, rng, pho, nuclides, d, dz, NL):
    while True:
        if g["E"]<E_CUT_G:# ④ 计数时 Dg_t 也要加这一份（审查文件第四节第 2 条）
            k=int(g["z"]/dz)
            k=min(k,NL-1)
            hist["Dg_c"][k]+=g["wt"]*g["E"]
            return
        sig=photon_sigma(pho,nuclides,g["E"])
        s=-math.log(1 - rng.random())/sig["tot"]
        z1=g["z"]+s*g["w"]
        if z1>d or z1<0:
            return 
        else:
            g["z"]=z1
        x=rng.random()*sig["tot"]
        if x<sig["pe"]:#光电
            k=int(g["z"]/dz)
            k=min(k,NL-1)
            hist["Dg_c"][k]+=g["wt"]*g["E"]
            return 
        elif x<sig["pe"]+sig["incoh"]:#康普顿⑤ KN 写好前先占位：光子直接结束，这份能量暂时没记
            mu=sample_kn_mu(g["E"],rng)
            E_new=g["E"]/(1+(g["E"]/ME_C2)*(1-mu))
            k=int(g["z"]/dz)
            k=min(k,NL-1)
            hist["Dg_c"][k]+=g["wt"]*(g["E"]-E_new)
            g["E"]=E_new
            phi=2*math.pi*rng.random()
            g["u"],g["v"],g["w"]=rotate_dir(g["u"],g["v"],g["w"],mu,phi)
            
        elif x<sig["pe"]+sig["incoh"]+sig["coh"]:#相干
            elem=pick_coh_elem(pho,nuclides,g["E"],sig["coh"],rng)
            mu=sample_coh_mu(pho[elem],g["E"],rng)
            phi=2*math.pi*rng.random()
            g["u"],g["v"],g["w"]=rotate_dir(g["u"],g["v"],g["w"],mu,phi)
        else:#电子对
            k=int(g["z"]/dz)
            k=min(k,NL-1)
            hist["Dg_c"][k]+=g["wt"]*(g["E"]-2*ME_C2)
            u,v,w=isotropic_direction(rng)
            dirc1={
                "z":g["z"],
                "u":u,
                "v":v,
                "w":w,
                "E":ME_C2,
                "wt":g["wt"]
            }
            hist["pstack"].append(dirc1)
            dirc2={
                    "z":g["z"],
                    "u":-u,
                    "v":-v,
                    "w":-w,
                    "E":ME_C2,
                    "wt":g["wt"]
                    }
            hist["pstack"].append(dirc2)
            return 

def run_photons(hist, rng, pho, nuclides, d, dz, NL):
    for site in hist["gsites"]:
        n=emit_photons(site,hist,rng)
    while hist["pstack"]:
        g=hist["pstack"].pop()
        track_photon(g,hist,rng,pho,nuclides,d,dz,NL)
    
def loglog0(x,y,e):
    assert x[0] <= e < x[-1], "loglog0：能量 e = %r eV 不在表 [%r, %r) 里 —— 光子掉到表下限以下还在查截面？1 keV 截断写在查截面之前了吗" % (e, x[0], x[-1])
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
            isf=g["coherent/integrated_scattering_factor"][()]
            assert isf.ndim == 2 and isf.shape[0] == 2, "integrated_scattering_factor 形状是 %r，不是 (2, N) —— 行 / 列的含义要重新核" % (isf.shape,)
            photon[elem]["coh_x2"]=isf[0,:]
            photon[elem]["coh_I"]=isf[1,:]
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

def pick_coh_elem(pho,nuclides,E,sig_coh,rng):
    r=rng.random()*sig_coh
    acc=0.0
    for nuc in nuclides:
        elem=NUC2ELEM[nuc["name"]]
        dens=nuc["dens"]
        acc+=loglog0(pho[elem]["E"],pho[elem]["coh"],E)*(1e-24)*dens
        if acc>r:
            return elem
    return elem

def sample_coh_mu(pe,E,rng):
    x2max=(E/HC)**2
    Imax=np.interp(x2max,pe["coh_x2"],pe["coh_I"])
    while True:
        r=rng.random()*Imax
        x2=np.interp(r,pe["coh_I"],pe["coh_x2"])
        mu=1-2*x2/x2max
        if rng.random()<0.5*(1+mu**2):
            return mu


def sample_kn_mu(E,rng):
    alpha=E/ME_C2
    while True:
        mu=2*rng.random()-1
        k=1/(1+alpha*(1-mu))
        f=k**2*(k+1/k-1+mu**2) 
        if rng.random()<f/2:
            return mu
    
def trap(y,h):
    return h*(y[0]/2+y[1:-1].sum()+y[-1]/2)

def kn_T_E(E):
    alpha=E/ME_C2
    mu=np.linspace(-1,1,2001)
    h=2/2000
    k=1/(1+alpha*(1-mu))
    f=k**2*(k+1/k-1+mu**2)
    return 1-(trap(f*k,h)/trap(f,h))
    
def emit_photons(site,hist,rng):
    xi=rng.random()
    n=int(xi+site["y"])
    for i in range(n):
        dic={}
        xi1=rng.random()
        if site["share"]>xi1:
            E_g=E_G_H
        else:
            E_g=E_G_B
        u,v,w=isotropic_direction(rng)
        dic["z"]=site["z"]
        dic["u"]=u
        dic["v"]=v
        dic["w"]=w
        dic["E"]=E_g
        dic["wt"]=site["wt"]
        hist["pstack"].append(dic)
    return n 

if __name__=="__main__":
    from water import kur_table,KEEP_FACTOR,build_water
    pho=load_photon()
    for elm in PHOTON_ELEMS:
        print(elm,len(pho[elm]["E"]),pho[elm]["E"][0],pho[elm]["E"][-1],len(pho[elm]["pair"]))
        print(pho[elm]["coh_x2"][0:3],pho[elm]["coh_x2"][-1])
        print(pho[elm]["coh_I"][0:3],pho[elm]["coh_I"][-1])
    edges,shares=kur_table()
    nuclides,order=build_water(KEEP_FACTOR*edges[-1])
    for E in [2.2233e6,1e6,0.4778e6]:
        sig=photon_sigma(pho,nuclides,E)
        print(E,sig["tot"],sig["pair"]/sig["tot"])
    for E in [1e4,1e6,2.2233e6]:
        print(kn_T_E(E))
