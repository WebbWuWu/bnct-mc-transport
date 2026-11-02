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
from check_thermal import check_T1,check_T1_selftest,check_T3E,check_T3E_selftest
from rng import lcg

SAB_PATH = os.path.join(r"E:\BNCT_GPU\lib80x_hdf5", "c_H_in_H2O.h5")   # 水里 H 的热散射库，约 629 MB
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
    return E_new

if __name__=="__main__":
    thermal=load_thermal(SAB_PATH,T_TAG)
    print(thermal["E_max"],len(thermal["xs_E"]),np.interp(0.0253,thermal["xs_E"], thermal["xs_sig"]))
    print(len(thermal["tables"]),thermal["tables"][48]["mu"].shape,thermal["tables"][48]["cdf"][-1],np.array_equal(thermal["xs_E"], thermal["E_in"]))
    # check_T1(thermal)
    # check_T1_selftest(thermal)
    # for E in [0.0265,0.0253,3e-6]:
    #     l,f=sample_thermal(thermal,E,None)
    #     print(l,f)
    rng=lcg(20261102)
    # for i in range(10):
    #     a,b,c=sample_thermal(thermal,0.0265,rng)
    #     print(a,b,c,(a-b),a/b)
    check_T3E(thermal, sample_thermal)
    check_T3E_selftest(thermal, sample_thermal)
