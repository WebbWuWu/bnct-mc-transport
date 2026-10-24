# ============================================================
#  rewrite_fgk_v2_1.py —— 默写 free_gas_kinematics 精简版（第 1 次，10/23）
#  2026-10-23 Claude 开文件，只写了这段注释；下面全部由你写
#  精简版 = 10/17 的抄写稿 whiteboard\copy_fgk_v2.py（三个分量装进一个长度 3 的列表，配 dot、norm 两个小函数）
# ============================================================
#  要写：
#    1. import
#    2. dot(a, b)、norm(a)
#    3. free_gas_kinematics(E, u, v, w, A, kT, x, mu_T, rng) → 返回 E', u', v', w'
#       步骤（只给步名，内容自己回忆；想不起来的空着写 # ?）：
#         0 算 y   1 中子速度   2a 抽辅助方向、去掉沿中子方向的部分（最多 100 轮）
#         2b 归一化成 e1   2c 靶核方向   2d 靶核速度   3 质心速度
#         4 CM 里的中子速度和速率   5 CM 里各向同性抽新方向   6 回 LAB → E' 和新方向
#  不写：assert、if __name__ == "__main__": 那一行和下面整块、print
#  不看：whiteboard\copy_fgk_v2.py、lesson07_freegas.py、讲解_free_gas_kinematics.md、whiteboard\rewrite_fgk1.py
#  写完：存盘，把文件名发给 Claude 对拍（同输入同种子，四个返回值 + 发生器状态逐位比）
# ============================================================
import math
from rng import lcg, isotropic_direction
import os
import h5py
import numpy as np
from xs_lookup import LIB,sigma_at,find_cell,load_nuclide

def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def norm(a):
    return math.sqrt(a[0]**2+a[1]**2+a[2]**2)

def free_gas_kinematics(E,u,v,w,A,kT,x,mu_T,rng):
    d=(u,v,w)
    y=math.sqrt(E*A/kT)
    vn=[y*d[i] for i in range(3)]#中子速度
    # 2a. 抽一个方向 p，减掉它沿 d 的部分
    for k in range(100):
        p=isotropic_direction(rng)
        al=dot(p,d)
        r=[p[i]-d[i]*al for i in range(3)]
        if norm(r)>=1e-6:
            break
    # 2b. r 归一化 → e1
    nr=norm(r)
    e1=[r[i]/nr for i in range(3)]
    # 2c. 靶核方向 = mu_T × d + sqrt(1 − mu_T²) × e1
    st=math.sqrt(1-mu_T**2)
    rb=[mu_T*d[i]+st*e1[i] for i in range(3)]
    # 2d. 靶核速度 = x × 靶核方向
    vb=[x*rb[i] for i in range(3)]
    # 3.  质心速度 = (中子速度 + A × 靶核速度) / (1 + A)
    vcm=[(vn[i]+A*vb[i])/(1+A) for i in range(3)]
    # 4.  CM 里中子速度 = 中子速度 − 质心速度，记下速率 s
    vc=[vn[i]-vcm[i] for i in range(3)]
    s=norm(vc)
    q= isotropic_direction(rng)
    vl=[vcm[i]+q[i]*s for i in range(3)]
    L=norm(vl)
    return (kT/A)*L**2,vl[0]/L,vl[1]/L,vl[2]/L
