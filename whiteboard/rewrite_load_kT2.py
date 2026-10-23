# ============================================================
#  rewrite_load_kT2.py —— 默写 load_kT + sigma_at_ext（第 2 次）
#  2026-10-19 Claude 开文件，只写了这段注释；下面全部由你写
# ============================================================
#  要写：
#    1. import（要哪些自己想；两个函数用不到的就不写）
#    2. load_kT(nuc, T)：打开一个核素的 HDF5 文件，读出温度标签 T（比如 "294K"）对应的 kT（eV），返回它
#    3. sigma_at_ext(E, XS, e)：查能量 e 处各道截面和总截面；e 低于能量表下限时按 1/v 往下外推
#       （返回哪几个值、外推的系数怎么来，自己回忆）
#  不写：assert、if __name__ == "__main__": 那一行和下面整块、print
#  不看：copy_load_kT.py、lesson07_freegas.py、讲解_load_kT.md、rewrite_load_kT1.py
#  水课在纸上手写过的：晚上照纸敲进来就行（照纸敲不算看原文）
#  写完：存盘，告诉 Claude 文件名 —— Claude 代对拍（容器没有 h5py：load_kT 只能逐行读，
#        sigma_at_ext 用假库对拍外推 / 不外推两个点逐位）
# ============================================================
import os
import math
import h5py
from xs_lookup import LIB,sigma_at

def load_kT(nuc,T):
    path=os.path.join(LIB,nuc+".h5")
    with h5py.File(path,"r")as f:
        root=f[nuc]
        kT=float(root["kTs"][T][()])
    return kT
def sigma_at_ext(E,XS,e):
    if e>=E[0]:
        XS,total=sigma_at(E,XS,e)
    else:
        XS1,total1=sigma_at(E,XS,E[0])
        XS=[0.0]*len(XS1)
        for i in range(len(XS1)):
            XS[i]=XS1[i]*math.sqrt(E[0]/e)
        total=total1*math.sqrt(E[0]/e)
    return XS,total
        