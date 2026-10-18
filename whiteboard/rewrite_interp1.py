# ============================================================
#  rewrite_interp1.py —— 默写 interp 三个（第 1 次）：find · interp_loglog · number_density
#  2026-10-15 Claude 开文件，只写了这段注释；下面是你 10/18 凭记忆默写的（实为第 2 次；10/15 纸稿没敲进来）
# ============================================================
#  要写：
#    1. import（要哪些自己想；函数本身用不到的就不写）
#    2. 三个常数 NA、M_N、RHO（名字和数值）
#    3. find（在能量表里找 e 落在哪一格，返回左端下标；两头越界怎么办）
#    4. interp_loglog（log-log 插值）
#    5. number_density（质量份额 + awr → 数密度）
#  不写：assert、if __name__ == "__main__": 那一行和下面整块、print
#  不看：copy_interp3.py、lesson05_tissue_scaffold.py、whiteboard_05_rewrite.py、
#        whiteboard\whiteboard_05_scaffold.py、讲解_interp三个.md
#  写完：存盘，把文件名发给 Claude 对拍（同一组输入跑原函数和你的函数，逐位比）
# ============================================================
import math
import numpy as np

NA   = 6.02214076e23
M_N  = 1.00866491588
RHO  = 1.06

def find(E,e):
    n=len(E)
    if e>E[-1]:
        return n-2
    if e<E[0]:
        return 0
    else:
        lo=0
        hi=n-1
        while True:
            mid=lo+(hi-lo)//2
            if e<E[mid]:
                hi=mid
            if e>E[mid]:
                lo=mid
            if e==E[mid]:
                return mid
            if hi-lo==1:
                return lo

def interp_loglog(E_tab,S_tab,e):
    k=int(np.searchsorted(E_tab,e,side="right"))-1
    k=min(max(k,0),len(E_tab)-2)
    p=np.log(S_tab[k+1]/S_tab[k])/np.log(E_tab[k+1]/E_tab[k])
    return S_tab[k]*(e/E_tab[k])**p

def number_density(w,awr):
    M=awr*M_N
    return w*RHO*NA/M
