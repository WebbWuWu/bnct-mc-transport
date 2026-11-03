# -*- coding: utf-8 -*-
# ============================================================
#  elastic_angle.py  --  M3：弹性散射（MT=2）真实角分布（质心系）
#  函数：load_elastic_angle（读库切表，11/3）· 以后加：按表抽 μ_cm、质心系 → 实验室系
#  文件头（import、常量）Claude 11/3 写；函数你写。依据：explore\lesson12_elastic_angle.py 的 show_angle
#  用法清单：项目文档\load_elastic_angle用法_2026-11-03.md
#  判据：Claude 写（写完第 3 步给）
#  运行：Windows 的 openmc (3.10.20) 环境，Ctrl+F5（只用 h5py + numpy）
# ============================================================

import os
import h5py
import numpy as np

LIB = r"E:\BNCT_GPU\lib80x_hdf5"     # 截面库目录（和 lesson12_elastic_angle.py 同一个）
NUCS = ("H1", "O16")                 # 水模体里要换真实角分布的两个核素

for _nuc in NUCS:
    assert os.path.exists(os.path.join(LIB, _nuc + ".h5")), "找不到 %s 的库文件" % _nuc

# ===== 第 1 步起从这里往下写 =====
def load_elastic_angle(lib_dir,nucs):
    angle={}
    for nuc in nucs:
        path=os.path.join(lib_dir,nuc+".h5")
        with h5py.File(path,"r")as f:
            cm=bool(f[nuc]["reactions"]["reaction_002"].attrs["center_of_mass"])
            ang=f[nuc]["reactions"]["reaction_002"]["product_0/distribution_0/angle"]
            ang_E=ang["energy"][:]
            mu_all=ang["mu"][:]
            mu_offs=ang["mu"].attrs["offsets"]
            n_pts=np.diff(np.append(mu_offs,mu_all.shape[1]))
            tables=[]
            for j in range(len(mu_offs)):
                b=mu_offs[j]
                e=b+n_pts[j]
                mu_v,mu_p,mu_c=mu_all[0:3,b:e]
                table1={
                    "mu":mu_v,
                    "pdf":mu_p,
                    "cdf":mu_c
                }
                tables.append(table1)
            angle[nuc]={"E":ang_E,"cm":cm,"tables":tables}
    return angle

if __name__ == "__main__":
    from check_elastic import check_A, check_A_selftest
    angle=load_elastic_angle(LIB, NUCS)
    check_A(angle)
    check_A_selftest(angle)
