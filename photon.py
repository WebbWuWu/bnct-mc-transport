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
from rng import lcg
from mcstat import mean_se
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
    n_mean=sig_cap/xs["Sigma_t"]          # 11/28 改名（原来叫 y）：每次碰撞平均生几个光子；"y" 让给位置的 y 坐标
    share=xs["gH_coe"]/sig_cap
    dic={
        "x":p.get("x",0.0),               # 11/28 M5 搬运 ②：多存 x、y（三维要用）；一维的 p 里没有 x、y，存 0.0
        "y":p.get("y",0.0),
        "z":p["z"],
        "wt":p["wt"],
        "n_mean":n_mean,
        "share":share,
    }
    hist["gsites"].append(dic)


def track_layer(arr,z0,z_end,w,s,val,dz,NL):
    assert 0<=z0<=NL*dz,"z0在板外,z0 = %g 不在 [0, %g] 里" % (z0, NL*dz)
    seg=0.0
    if w==0:
        k0=min(int(z0/dz),NL-1)
        arr[k0]+=val*s
        seg=s
        expect=s
    else:
        k0=min(int(z0/dz),NL-1)
        k1=min(int(z_end/dz),NL-1)
        cur=z0
        k=k0
        guard=0
        while True:
            guard+=1
            assert guard<=NL+2,\
            "拆层循环未退出"
            if w>0.0:
                k_o=(k+1)*dz
                z_step=min(z_end,k_o)
            else:
                k_o=k*dz
                z_step=max(z_end,k_o)
            assert 0<=k<NL,\
                "k越界"
            dz_step=abs(z_step-cur)
            arr[k]+=val*dz_step/abs(w)
            seg+=dz_step/abs(w)
            cur=z_step
            if k==k1:
                break
            elif w>0.0:
                k+=1
            else:
                k-=1
        expect=abs(z_end-z0)/abs(w)
    assert math.isclose(seg, expect, rel_tol=1e-9, abs_tol=1e-12), \
    "段记账不守恒: seg=%.15f expect=%.15f z0=%.6f z_end=%.6f w=%.9f" \
    % (seg, expect, z0, z_end, w)
    return seg

def track_photon(g, hist, rng, pho, nuclides, d, dz, NL):
    while True:
        if g["E"]<E_CUT_G:# ④ 计数时 Dg_t 也要加这一份（审查文件第四节第 2 条）
            k=int(g["z"]/dz)
            k=min(k,NL-1)
            hist["Dg_c"][k]+=g["wt"]*g["E"]
            hist["Dg_t"][k]+=g["wt"]*g["E"]
            return
        sig=photon_sigma(pho,nuclides,g["E"])
        s=-math.log(1 - rng.random())/sig["tot"]
        z1=g["z"]+s*g["w"]
        if z1>d:
            z_end=d
        elif z1<0:
            z_end=0.0
        else:
            z_end=z1
        mt=mu_tr(sig,g["E"],pho)
        track_layer(hist["Dg_t"],g["z"],z_end,g["w"],s,g["wt"]*mt*g["E"],dz,NL)
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

def run_beam(N,seed,E0,pho,nuclides,d,NL):
    dz=d/NL
    rng=lcg(seed)
    t_sum=[0.0]*NL
    t_sq=[0.0]*NL
    c_sum=[0.0]*NL
    c_sq=[0.0]*NL
    for i in range(N):
        hist={
            "Dg_c":[0.0]*NL,
            "Dg_t":[0.0]*NL,
            "pstack":[],
        }
        g={
            "z":0.0,
            "u":0.0,
            "v":0.0,
            "w":1.0,
            "E":E0,
            "wt":1.0
        }
        hist["pstack"].append(g)
        while hist["pstack"]:
            g=hist["pstack"].pop()
            track_photon(g,hist,rng,pho,nuclides,d,dz,NL)
        for k in range(NL):
            t_sum[k]+=hist["Dg_t"][k]
            t_sq[k]+=hist["Dg_t"][k]**2
            c_sum[k]+=hist["Dg_c"][k]
            c_sq[k]+=hist["Dg_c"][k]**2
    t_mean=[]
    t_se=[]
    c_mean=[]
    c_se=[]
    for k in range(NL):
        mean,_,se=mean_se(t_sum[k],t_sq[k],N)
        t_mean.append(mean)
        t_se.append(se)
        mean1,_,se1=mean_se(c_sum[k],c_sq[k],N)
        c_mean.append(mean1)
        c_se.append(se1)
    return {
        "t_mean":t_mean,
        "t_se":t_se,
        "c_mean":c_mean,
        "c_se":c_se,
        "N":N,
    }

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
    photon["KN"]={}
    photon["KN"]["E"]=np.geomspace(E_CUT_G,2.3e6,200)
    photon["KN"]["TE"]=np.zeros(len(photon["KN"]["E"]))
    for i in range(len(photon["KN"]["E"])):
        photon["KN"]["TE"][i]=kn_T_E(photon["KN"]["E"][i])
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

def kn_TE_tab(pho,E):
    Eg=pho["KN"]["E"]
    TE=pho["KN"]["TE"]
    assert Eg[0] <= E<=Eg[-1],\
        "E越界,E=%e,Eg[0]=%e,Eg[-1]=%e"%(E,Eg[0],Eg[-1])
    return np.interp(E,Eg,TE)

def mu_tr(sig,E,pho):
    return sig["pe"]+sig["incoh"]*kn_TE_tab(pho,E)+sig["pair"]*(1-2*ME_C2/E)
    
def emit_photons(site,hist,rng):
    xi=rng.random()
    n=int(xi+site["n_mean"])        # 11/28：键名 "y" 改成 "n_mean"（M5 搬运 ②）
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
        # [Claude 判据] μ_tr vs NIST 水的 μ_en/ρ（另一条独立的路）。
    # 我们的水 ρ = 0.998，所以参考值 = NIST × 0.998；硼只有 15 ppm，可以忽略。
    # 门槛（看数之前定）：相对差 < 3%（差在纯 KN 角分布、μ_tr 和 μ_en 差一个辐射份额 g ≈ 0.1%）
    NIST_MUEN = {1e5: 0.02546, 1e6: 0.03103, 2e6: 0.02608}   # cm²/g，NIST 水
    for E in [1e5, 1e6, 2e6]:
        sig = photon_sigma(pho, nuclides, E)
        mt = mu_tr(sig, E, pho)
        ref = NIST_MUEN[E] * 0.998
        print("mu_tr E=%g eV: 程序 %.5f /cm, 参考 %.5f /cm, 相对差 %+.2f%%" % (E, mt, ref, 100 * (mt / ref - 1)))
        assert abs(mt / ref - 1) < 0.03, "mu_tr 在 E=%g eV 和 NIST 差超过 3%%：程序 %.5f，参考 %.5f" % (E, mt, ref)
        