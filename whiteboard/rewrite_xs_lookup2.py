# ============================================================
#  rewrite_xs_lookup2.py —— 默写 xs_lookup 四个函数（第 2 次，10/23）
#  2026-10-23 Claude 开文件，只写了这段注释；下面全部由你写
#  放在 whiteboard\ 下：import 自己的模块时不用 —— 这个文件只用标准库 + numpy + h5py
# ============================================================
#  要写：
#    1. import + 两个常数：LIB（截面库目录）、BARN
#    2. load_nuclide(nuc, T)：读一个核素 → 返回 E, mts, XS, awr, Q（跳过 redundant 的道；按阈值下标补 0）
#    3. xs_at(E, s, e, lo)：一个道在 e 处线性插值
#    4. find_cell(E, e)：e 落在第几个格子（含 e == E[-1] 的特殊情况）
#    5. sigma_at(E, XS, e)：所有道在 e 处的截面列表 + 总截面
#  不写：assert、if __name__ == "__main__": 那一行和下面整块、print
#  不计分：NUC、T、E_TEST 这几个只给 __main__ 用的常数（写不写都行）
#  不看：xs_lookup.py、copy_xs_lookup.py、讲解_xs_lookup.md、lesson11_mt301.py
#  写完：存盘，把文件名发给 Claude 对拍（假库查值逐位比；load_nuclide 我这边没有 h5py，只做逐行读）
# ============================================================
import os
import h5py
import numpy as np

LIB=r"E:\BNCT_GPU\lib80x_hdf5"
NUC="B10"
T="294K"
E_TEST=0.0253
BARN=1.0e-24

def load_nuclide(nuc,T):
    path=os.path.join(LIB,nuc+".h5")
    with h5py.File(path,"r")as f:
        root=f[nuc]
        E=root["energy"][T][:]
        awr=float(root.attrs["atomic_weight_ratio"])
        names=list(root["reactions"].keys())
        mts=[]
        XS=[]
        Q=[]
        for name in sorted(names):
            redundant=root["reactions"][name].attrs["redundant"]
            if redundant==1:
                continue
            xs_ds=root["reactions"][name][T]["xs"]
            xs=xs_ds[:]
            thr=int(xs_ds.attrs["threshold_idx"])
            s=np.zeros(len(E))
            s[thr:]=xs
            XS.append(s)
            mts.append(int(root["reactions"][name].attrs["mt"]))
            Q.append(float(root["reactions"][name].attrs["Q_value"]))
    return E,mts,XS,awr,Q

def xs_at(E,s,e,lo):
    i=lo+1
    frac=(e-E[lo])/(E[i]-E[lo])
    return s[lo]+(s[i]-s[lo])*frac

def find_cell(E,e):
    if e==E[-1]:
        lo=len(E)-2
    else:
        lo=np.searchsorted(E,e,side="right")-1
    return lo 

def sigma_at(E,XS,e):
    lo=find_cell(E,e)
    sig=[]
    for i in range(len(XS)):
        sig.append(xs_at(E,XS[i],e,lo))
    total=sum(sig)
    return sig,total
