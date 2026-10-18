# ============================================================
#  rewrite_total_xs1.py —— 默写 total_xs（第 1 次）
#  2026-10-15 Claude 开文件，只写了这段注释；下面是你 10/18 凭记忆默写的（实为第 2 次；10/15 纸稿没敲进来）
# ============================================================
#  要写：
#    1. import（要哪些自己想；函数本身用不到的就不写）
#    2. 三个常数 LIB（库路径）、T（温度标签）、E_TEST（函数里要用）—— 路径字符串照抄、不计分
#    3. total_xs（打开一个核素的 HDF5 文件，把该加的道加成总截面；返回哪几个值自己回忆）
#  不写：assert、if __name__ == "__main__": 那一行和下面整块、print
#  不看：copy_total_xs.py、xslib.py、讲解_total_xs.md
#  写完：存盘，把文件名发给 Claude 对拍（在你电脑上跑原函数和你的函数，逐位比）
# ============================================================
import os
import h5py
import numpy as np

LIB    = r"E:\BNCT_GPU\lib80x_hdf5"
T      = "294K"
E_TEST = 0.0253

def total_xs(nuc,T):
    path=os.path.join(LIB,nuc+".h5")
    with h5py.File(path,"r")as f:
        root=f[nuc]
        E=root["energy"][T][:]
        awr=float(root.attrs["atomic_weight_ratio"])
        names=list(root["reactions"].keys())
        total=np.zeros(len(E))
        at_i0={}
        i0=int(np.searchsorted(E,E_TEST,side="right"))-1
        for name in sorted(names):
            rx=root["reactions"][name]
            mt=rx.attrs["mt"]
            redundant=rx.attrs["redundant"]
            xs_ds=rx[T]["xs"]
            thr=int(xs_ds.attrs["threshold_idx"])
            xs=xs_ds[:]
            at_i0[int(mt)]=0.0 if i0<thr else xs[i0-thr]
            if redundant:
                continue
            total[thr:]+=xs
    return E,total,at_i0,i0,awr
