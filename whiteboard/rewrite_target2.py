# ============================================================
#  rewrite_target2.py —— 默写 sample_target（第 2 次）
#  2026-10-11 Claude 开文件，只写了这段注释；下面全部由你写
# ============================================================
#  要写：
#    1. import（要哪些自己想）
#    2. path_A_prob
#    3. sample_target
#  不写：assert、if __name__ == "__main__": 那一行和下面整块、print
#  不看：lesson07_freegas.py、rewrite_target1.py、rewrite_sample.py、讲义
#  写完：存盘，把文件名发给 Claude 对拍（同种子跑原函数和你的函数，逐位比）
# ============================================================
import math
import numpy as np
from rng import lcg

def path_A_prob(y):
    return y*math.sqrt(math.pi)/(y*math.sqrt(math.pi)+2)

def sample_target(y,rng):
    for k in range(1000):
        xi2=rng.random()
        xi3=rng.random()
        xi4=rng.random()
        if rng.random()<path_A_prob(y):
            x2=-np.log(1-xi2)-np.log(1-xi3)*math.cos((1-xi4)/2*math.pi)**2
        else:
            x2=-np.log((1-xi2)*(1-xi3))
        x=math.sqrt(x2)
        mu_T=2*rng.random()-1
        v_rel=math.sqrt(max((x2+y**2-2*x*y*mu_T),0.0))
        pac=v_rel/(x+y)
        if rng.random()<pac:
            return x,mu_T,k+1
