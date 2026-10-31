# -*- coding: utf-8 -*-
# ============================================================
#  lesson12_sab_explore.py  --  第 12 课：c_H_in_H2O.h5 结构探索（只看「目录」，不读大数组）
#  文件头（import、路径常量、存在性 assert）Claude 10/31 写；第 2 步起你写。
#  运行：Windows 的 openmc (3.10.20) 环境，Ctrl+F5。只用 h5py + numpy，不用 OpenMC，不进 WSL。
#  用法：项目文档\c_H_in_H2O探索用法_2026-10-30.md
# ============================================================

import os
import h5py
import numpy as np

LIB = r"E:\BNCT_GPU\lib80x_hdf5"
SAB_PATH = os.path.join(LIB, "c_H_in_H2O.h5")    # 水里 H 的热散射库，约 629 MB
H1_PATH = os.path.join(LIB, "H1.h5")              # 一致性检查用：H1 自己的 294K kT
K_B = 8.617333e-5                                 # 玻尔兹曼常数，eV/K（kT ÷ K_B = 温度 K）
T_TAG = "294K"                                    # 要看的温度档名字

assert os.path.exists(SAB_PATH), "找不到热散射库：%s" % SAB_PATH
assert os.path.exists(H1_PATH), "找不到 H1 库：%s" % H1_PATH

# ===== 第 2 步起从这里往下写 =====
def show_item(path,obj):
    if isinstance(obj,h5py.Dataset):
        print(path,"数据集",obj.shape,obj.dtype)
    else:
        print(path,"组")

with h5py.File(SAB_PATH,"r")as f:
    names=list(f.keys())
    for name in names:
        print("%s"%name)
    root=f["c_H_in_H2O"]
    for key,val in root.attrs.items():
        print(key,val)
    kts=root["kTs"]
    for kt in list(kts.keys()):
        print(kt)
    kT_sab=float(kts[T_TAG][()])
    print(kT_sab,kT_sab/K_B)
    tgrp=root[T_TAG]
    tgrp.visititems(show_item)
    xs_ds=tgrp["inelastic/xs"]
    print(xs_ds.shape)
    for key,val in xs_ds.attrs.items():
        print(key,val)
    xs_arr=xs_ds[:]
    E_grid=xs_arr[0]
    sig_grid=xs_arr[1]
    print("E_grid第一个数=%e,E_grid最后一个数=%e"%(E_grid[0],E_grid[-1]))
    print("能量=%e，截面=%e"%(0.0253,np.interp(0.0253,E_grid, sig_grid)))
    print("能量=%e，截面=%e"%(1,np.interp(1,E_grid, sig_grid)))
    print("能量=%e，截面=%e"%(10,np.interp(10,E_grid, sig_grid)))
    dist=tgrp["inelastic/distribution"]
    E_in=dist["energy"][:]
    j0=np.argmin(np.abs(E_in-0.0253))
    print("j0=%d,E_in[j0]=%e,E_in[j0-1]=%e,E_in[j0+1]=%e"%(j0,E_in[j0],E_in[j0-1],E_in[j0+1]))
    eout=dist["energy_out"][:]
    offs=dist["energy_out"].attrs["offsets"]
    i_beg=offs[j0]
    i_end=offs[j0+1]
    E_out,pdf_out,cdf_out=eout[0:3,i_beg:i_end]
    print("第 %d 张表：入射能量 %e eV，点数 %d" % (j0, E_in[j0], i_end - i_beg))
    print("E_out 第一个 = %e，最后一个 = %e" % (E_out[0], E_out[-1]))
    print("cdf_out 第一个 = %e，最后一个 = %e" % (cdf_out[0], cdf_out[-1]))
    print("pdf_out 对 E_out 的积分 = %e" % np.trapezoid(pdf_out, E_out))
    mu_interp,mu_off=eout[3:5,i_beg:i_end]
    print("mu_interp_max=%e,mu_interp_min=%e"%(max(mu_interp),min(mu_interp)))
    print("mu_off_max=%e,mu_off_min=%e"%(max(mu_off),min(mu_off)))
    mu_all=dist["mu"][:]
    m_beg=int(mu_off[300])
    m_end=int(mu_off[301])
    mu_val,mu_pdf,mu_cdf=mu_all[0:3,m_beg:m_end]
    print("第 300 个出射点：E_out = %e eV，μ 表点数 %d" % (E_out[300], m_end - m_beg))
    print("mu_val =", mu_val)
    print("mu_pdf =", mu_pdf)
    print("mu_cdf =", mu_cdf)
    for E_tgt in (0.0253,1.0,10.0):
        j=np.argmin(np.abs(E_in-E_tgt))
        b=offs[j]
        if j+1==len(offs):
            e=eout.shape[1]
        else:
            e=offs[j+1]    
        Eo,po,co=eout[0:3,b:e]
        f_up=1-np.interp(E_in[j],Eo,co)
        E_mean=np.trapezoid(Eo*po,Eo)
        print("入射 %e eV（第 %d 张）：E_out 范围 %e ~ %e，上散射比例 %.3f，平均出射能量/入射 = %.3f"
            % (E_in[j], j, Eo[0], Eo[-1], f_up, E_mean / E_in[j]))
    print("--- dist 的属性")
    for key,val in dist.attrs.items():
        print(key,val)
    print("--- energy 的属性")
    for key,val in dist["energy"].attrs.items():
        print(key,val)
    print("--- energy_out 的属性")
    for key,val in dist["energy_out"].attrs.items():
        print(key,val)
    print("--- mu 的属性")
    for key,val in dist["mu"].attrs.items():
        print(key,val)
    
with h5py.File(H1_PATH,"r")as g:
    kT_h1=float(g["H1"]["kTs"][T_TAG][()])
    print("kT_sab=%e,kT_h1=%e,(kT_sab − kT_h1) ÷ kT_h1=%e"%(kT_sab,kT_h1,(kT_sab-kT_h1)/kT_h1))
