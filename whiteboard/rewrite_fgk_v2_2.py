# ============================================================
#  rewrite_fgk_v2_2.py —— 默写 free_gas_kinematics 精简版（第 1 次，10/25）
#  2026-10-25 Claude 开文件，只写了这段注释；下面全部由你写
#  （编号说明：rewrite_fgk_v2_1.py 是 10/24 你照抄写稿敲的那一份，今天这份是第一次不看原文默）
# ============================================================
#  要写：
#    1. import（函数本身用得到的才写：开方要哪个库？各向同性方向从哪个文件来？）
#    2. dot(a, b)、norm(a)
#    3. free_gas_kinematics(E, u, v, w, A, kT, x, mu_T, rng) → 返回 E', u', v', w'
#       步骤（只给步名，内容自己回忆；想不起来的空着写 # ?）：
#         0 中子方向 d、算 y   1 中子速度
#         2a 抽辅助方向、去掉沿中子方向的部分（最多 100 轮，够长就停）
#         2b 归一化成 e1   2c 靶核方向   2d 靶核速度   3 质心速度
#         4 CM 里的中子速度和速率   5+6 CM 里各向同性抽新方向，回 LAB → E' 和新方向
#  不写：assert、if __name__ == "__main__": 那一行和下面整块、print
#        2a 那个 for 循环的 else 里只有一条 assert —— 所以 else 整块（连 else: 那一行）都不写
#  不看：whiteboard\copy_fgk_v2.py、whiteboard\rewrite_fgk_v2_1.py、whiteboard\rewrite_fgk1.py、
#        lesson07_freegas.py、项目文档\讲解_free_gas_kinematics.md
#  写完：存盘，把文件名发给 Claude 对拍（同输入同种子，四个返回值 + 发生器状态逐位比）
# ============================================================
import math
from rng import lcg, isotropic_direction

def dot(a,b):
    s=0.0
    for i in range(len(a)):
        s+=a[i]*b[i]
    return s

def norm(a):
    s=0.0
    for i in range(len(a)):
        s+=a[i]**2
    return math.sqrt(s)

def free_gas_kinematics(E,u,v,w,A,kT,x,mu_T,rng):
    d=(u,v,w)
    y=math.sqrt(E*A/kT)
    vn=[y*d[i] for i in range(len(d))]
    for k in range(100):
        p=isotropic_direction(rng)
        e=[p[i]-dot(p,d)*d[i] for i in range(len(p))]
        if norm(e)>=1e-6:
            break
    e1=[e[i]/norm(e)for i in range(len(e))]
    rb=[mu_T*d[i]+math.sqrt(1-mu_T**2)*e1[i]for i in range(len(e))]
    vb=[x*rb[i] for i in range(len(e))]
    vcm=[(vn[i]+A*vb[i])/(1+A)for i in range(len(e))]
    vcmn=[vn[i]-vcm[i]for i in range(len(e))]
    z=isotropic_direction(rng)
    ycmn=norm(vcmn)
    vl=[z[i]*ycmn+vcm[i]for i in range(len(e))]
    l=norm(vl)
    return (kT/A)*l**2,vl[0]/l,vl[1]/l,vl[2]/l
