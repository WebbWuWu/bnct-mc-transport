# -*- coding: utf-8 -*-
# ============================================================
#  thermal.py  --  M3：水里 H 的热散射 S(α,β)（c_H_in_H2O）
#  三个函数：load_thermal（读库切表）· sample_thermal（抽 E′、μ）· rotate_dir（转方向）
#  文件头（import、常量）Claude 11/1 写；函数你写。框架：项目文档\M3框架_S_alpha_beta_2026-10-31.md
#  判据 T1：check_thermal.py（Claude 写）
#  运行：Windows 的 openmc (3.10.20) 环境，Ctrl+F5（只用 h5py + numpy）
# ============================================================

import os,math
import h5py
import numpy as np
from rng import lcg
from lib_path import LIB_DIR

_LIB_DIR = LIB_DIR                                                   # 11/10：库路径统一写在 lib_path.py（值没变）
SAB_PATH = os.path.join(_LIB_DIR, "c_H_in_H2O.h5")                    # 水里 H 的热散射库，约 629 MB
T_TAG = "294K"                                                         # 只用 294K 这一档（M3 框架 ⑤）
N_MU=20

assert os.path.exists(SAB_PATH), "找不到热散射库：%s" % SAB_PATH

# ===== 第 1 步起从这里往下写 =====
def load_thermal(sab_path,T_tag):
    thermal={}
    with h5py.File(sab_path,"r")as f:
        root=f["c_H_in_H2O"]
        thermal["E_max"]=float(root.attrs["energy_max"])
        xs=root[T_tag]["inelastic/xs"][:]
        thermal["xs_E"]=xs[0,:]
        thermal["xs_sig"]=xs[1,:]
        dist=root[T_tag]["inelastic/distribution"]
        thermal["E_in"]=dist["energy"][:]
        eout=dist["energy_out"][:]
        offs=dist["energy_out"].attrs["offsets"]
        mu_all=dist["mu"][:]
        n_out=np.diff(np.append(offs,eout.shape[1]))
        thermal["tables"]=[]
        for l in range(len(offs)):
            beg=offs[l]
            end=beg+n_out[l]
            E_out,pdf,cdf=eout[0:3,beg:end]
            mu_beg=eout[4,beg:end].astype(int)
            mu_rows=[]
            for k in range(len(mu_beg)):
                mu_rows.append(mu_all[0,mu_beg[k]:mu_beg[k] + N_MU])
            mu_tab = np.array(mu_rows)
            dic={
                "E_out":E_out,
                "pdf":pdf,
                "cdf":cdf,
                "mu":mu_tab,
            }
            thermal["tables"].append(dic)
    return thermal

def sample_thermal(thermal,E,rng):
    E_in=thermal["E_in"]
    i=np.searchsorted(E_in,E,side="right")-1
    if i<0:
        i=0
    if i>len(E_in)-2:
        i=len(E_in)-2
    f=(E-E_in[i])/(E_in[i+1]-E_in[i])
    if f>0.5:
        l=i+1
    else:
        l=i
    E_l=E_in[l]
    tab=thermal["tables"][l]
    cdf=tab["cdf"]
    xi1=rng.random()
    j=np.searchsorted(cdf,xi1,side="right")-1
    j=min(max(0,j),len(cdf)-2)
    E_out=tab["E_out"]
    pdf=tab["pdf"]
    E_j=E_out[j]
    E_j1=E_out[j+1]
    p_j=pdf[j]
    p_j1=pdf[j+1]
    c_j=cdf[j]
    c_j1=cdf[j+1]
    g=(xi1-c_j)/(c_j1-c_j)
    s=(p_j1-p_j)/(E_j1-E_j)
    if s==0.0:
        d=(xi1-c_j)/p_j
    else:
        d=(math.sqrt(max(0.0,(p_j**2 + 2*s*(xi1-c_j))))-p_j)/s
    E_new=E_j+d

    if E_new<E_l/2:
        E_new=E_new*(2*E/E_l-1)
    else:
        E_new+=E-E_l
    xi2=rng.random()
    k=int(xi2*N_MU)
    mu_tab=tab["mu"]
    mu=mu_tab[j,k]+g*(mu_tab[j+1,k]-mu_tab[j,k])
    if k==0:
        mu_left=-1-(mu+1)
    else:
        mu_left=mu_tab[j,k-1]+g*(mu_tab[j+1,k-1]-mu_tab[j,k-1])
    if k==N_MU-1:
        mu_right=1+1-mu
    else:
        mu_right=mu_tab[j,k+1]+g*(mu_tab[j+1,k+1]-mu_tab[j,k+1])
    xi3=rng.random()
    width=min((mu-mu_left),(mu_right-mu))
    mu+=width*(xi3-0.5)
    return E_new,mu

def rotate_dir(u,v,w,mu,phi):
    cos_phi=math.cos(phi)
    sin_phi=math.sin(phi)
    a=math.sqrt(max(0.0,1-mu**2))
    b=math.sqrt(max(0.0,1-w**2))
    if b>1e-10:
        u_new=u*mu+a*(u*w*cos_phi-v*sin_phi)/b
        v_new=mu*v+a*(v*w*cos_phi+u*sin_phi)/b
        w_new=mu*w-a*b*cos_phi
        nrm=math.sqrt(u_new**2+v_new**2+w_new**2)
        return u_new/nrm,v_new/nrm,w_new/nrm
    else:
        b=math.sqrt(1-v**2)
        u_new=mu*u+a*(u*v*cos_phi+w*sin_phi)/b
        v_new=mu*v-a*b*cos_phi
        w_new=mu*w+a*(v*w*cos_phi-u*sin_phi)/b
        nrm=math.sqrt(u_new**2+v_new**2+w_new**2)
        return u_new/nrm,v_new/nrm,w_new/nrm

if __name__=="__main__":
    import sys                                                         # 11/10：判据挪进 checks\，只有直接跑本文件时才 import 它（以前写在文件开头）
    sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),"checks"))
    from check_thermal import check_T1,check_T1_selftest,check_T3E,check_T3E_selftest,check_T3,check_T3_selftest,check_T2,check_T2_selftest
    thermal=load_thermal(SAB_PATH,T_TAG)
    print(thermal["E_max"],len(thermal["xs_E"]),np.interp(0.0253,thermal["xs_E"], thermal["xs_sig"]))
    print(len(thermal["tables"]),thermal["tables"][48]["mu"].shape,thermal["tables"][48]["cdf"][-1],np.array_equal(thermal["xs_E"], thermal["E_in"]))
    # check_T1(thermal)
    # check_T1_selftest(thermal)
    # for E in [0.0265,0.0253,3e-6]:
    #     l,f=sample_thermal(thermal,E,None)
    #     print(l,f)
    rng=lcg(20261102)
    # check_T3(thermal, sample_thermal)
    # check_T3_selftest(thermal, sample_thermal)
    # check_T3E(thermal, sample_thermal)
    # check_T3E_selftest(thermal, sample_thermal)
    check_T2(rotate_dir)
    check_T2_selftest()
