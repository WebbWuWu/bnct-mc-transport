# -*- coding: utf-8 -*-
# ============================================================
#  elastic_angle.py  --  M3：弹性散射（MT=2）真实角分布（质心系）
#  函数：load_elastic_angle（读库切表，11/3）· 以后加：按表抽 μ_cm、质心系 → 实验室系
#  文件头（import、常量）Claude 11/3 写；函数你写。依据：explore\lesson12_elastic_angle.py 的 show_angle
#  用法清单：项目文档\load_elastic_angle用法_2026-11-03.md
#  判据：Claude 写（写完第 3 步给）
#  运行：Windows 的 openmc (3.10.20) 环境，Ctrl+F5（只用 h5py + numpy）
# ============================================================

import os,math
import h5py
import numpy as np

from lib_path import LIB_DIR
LIB = LIB_DIR                        # 截面库目录；11/10：统一写在 lib_path.py（值没变）
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

def sample_elastic_mu(ang,E,rng):
    E_g=ang["E"]
    i=np.searchsorted(E_g,E,side="right")-1
    i=min(max(0,i),len(E_g)-2)
    f=((E-E_g[i])/(E_g[i+1]-E_g[i]))
    f=min(max(0.0,f),1.0)
    xi1=rng.random()
    if f>xi1:
        l=i+1
    else:
        l=i
    tab=ang["tables"][l]
    mu_v=tab["mu"]
    pdf=tab["pdf"]
    cdf=tab["cdf"]
    xi2=rng.random()
    j=np.searchsorted(cdf,xi2,side="right")-1
    j=min(max(0,j),len(cdf)-2)
    mu_j=mu_v[j]
    mu_j1=mu_v[j+1]
    p_j=pdf[j]
    p_j1=pdf[j+1]
    c_j=cdf[j]
    s=(p_j1-p_j)/(mu_j1-mu_j)
    if s==0.0:
        d=(xi2-c_j)/p_j
    else:
        d=(math.sqrt(max(0.0,(p_j**2+2*s*(xi2-c_j))))-p_j)/s
    mu=mu_j+d
    mu=min(max(-1.0,mu),1.0)
    return mu
        
if __name__ == "__main__":
    import sys                                                         # 11/10：判据挪进了 checks\
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "checks"))
    from check_elastic import check_A, check_A_selftest
    angle=load_elastic_angle(LIB, NUCS)
    from rng import lcg
    # check_A(angle)
    # check_A_selftest(angle)
    rng=lcg(20261104)
    # ang={
    #     "E":[1.0,2.0],
    #     "tables":[
    #         {"mu":[-0.6, 0.2, 0.9],
    #         "pdf":[0.35, 0.80, 0.55],
    #         "cdf":[0, 0.46, 0.9325]
    #         },
    #         {"mu":[-0.6, 0.2, 0.9],
    #         "pdf":[0.35, 0.80, 0.55],
    #         "cdf":[0, 0.46, 0.9325]
    #         }
    #     ]
    # }
    # for i in range(5):
    #     a=sample_elastic_mu(ang,1.5,rng)
    #     print(a)
    from check_elastic import check_B, check_B_selftest
    check_B(angle, sample_elastic_mu)
    check_B_selftest(angle)
