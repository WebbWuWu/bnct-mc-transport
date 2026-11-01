# -*- coding: utf-8 -*-
# ============================================================
#  lesson12_elastic_angle.py  --  M3 第二条线：弹性散射（MT=2）质心系角分布，读库探索
#  文件头（import、路径常量、存在性 assert）Claude 11/1 写；第 1 步起你写。
#  运行：Windows 的 openmc (3.10.20) 环境，Ctrl+F5。只用 h5py + numpy，不用 OpenMC，不进 WSL。
#  今天只「读」和「看数」，不抽样、不改 run_slab。
# ============================================================

import os
import h5py
import numpy as np

LIB = r"E:\BNCT_GPU\lib80x_hdf5"
NUCS = ("H1", "O16")                       # 水模体里要换角分布的两个核素（B10 以后再加）
E_PICKS = (0.0253, 1.0e4, 5.0e5, 2.0e6)    # 要看的四个入射能量，eV：热点 · 10 keV · 0.5 MeV · 2 MeV（KUR 源上限）

for nuc in NUCS:
    assert os.path.exists(os.path.join(LIB, nuc + ".h5")), "找不到 %s 的库文件" % nuc

# ===== 第 1 步起从这里往下写 =====
def show_item(path,obj):
    if isinstance(obj,h5py.Dataset):
        print(path,"数据集",obj.shape,obj.dtype)
    else:
        print(path,"组")

def show_angle(nuc):
    p1=os.path.join(LIB,nuc+".h5") 
    with h5py.File(p1,"r")as f:
        rx=f[nuc]["reactions"]["reaction_002"]
        print(nuc,rx.attrs["center_of_mass"])
        ang=rx["product_0/distribution_0/angle"]
        ang_E=ang["energy"][:]
        mu_all=ang["mu"][:]
        mu_offs=ang["mu"].attrs["offsets"]
        mu_interp=ang["mu"].attrs["interpolation"]
        n_pts=np.diff(np.append(mu_offs,mu_all.shape[1]))
        print(len(ang_E),ang_E[0],ang_E[-1],len(mu_offs),mu_offs[0:5],np.unique(mu_interp),np.unique(n_pts),sum(n_pts))
        # [Claude 判据 A1] 抓：offsets 和入射能量对不上（张数错）
    assert len(mu_offs) == len(ang_E), "%s：offsets %d 个，入射能量 %d 个" % (nuc, len(mu_offs), len(ang_E))
    # [Claude 判据 A4] 抓：插值方式不是线性-线性 —— mubar 的辛普森公式只在 pdf 每段是直线时精确
    assert np.all(mu_interp == 2), "%s：有表的插值方式不是 2（线性-线性）：%s" % (nuc, np.unique(mu_interp))
    for jj in range(len(ang_E)):
        mv, mp, mc = mu_all[0:3, mu_offs[jj]:mu_offs[jj] + n_pts[jj]]
        # [Claude 判据 A2] 抓：切片起止错一格（会把下一张表的 −1 切进来，μ 就不再从 −1 单调升到 +1）
        assert mv[0] == -1.0 and mv[-1] == 1.0 and np.all(np.diff(mv) > 0), \
            "%s 第 %d 张：μ 不是从 −1 单调升到 +1（首 %g、尾 %g）" % (nuc, jj, mv[0], mv[-1])
        # [Claude 判据 A3] 抓：pdf / cdf 行拿错（两条路：对 pdf 做梯形积分 vs 库里 cdf 的增量）
        assert abs(np.trapezoid(mp, mv) - (mc[-1] - mc[0])) < 1e-5, \
            "%s 第 %d 张：∫pdf = %.8f，cdf 增量 = %.8f" % (nuc, jj, np.trapezoid(mp, mv), mc[-1] - mc[0])
    print("%s 判据 A1–A4 过：%d 张表" % (nuc, len(ang_E)))
    for E_tgt in E_PICKS:
        j=np.argmin(np.abs(E_tgt-ang_E)) 
        b=mu_offs[j]
        e=b+n_pts[j]
        mu_v,mu_p,mu_c=mu_all[0:3,b:e]
        norm=np.trapezoid(mu_p,mu_v)
        mubar=np.sum((mu_v[1:]-mu_v[:-1])/6*((2*mu_v[:-1]+mu_v[1:])*mu_p[:-1]+(mu_v[:-1]+2*mu_v[1:])*mu_p[1:]))
        print(j,ang_E[j],n_pts[j],mu_v,mu_p,mu_c[-1],norm,mubar) 

for nuc in NUCS:
    show_angle(nuc)
