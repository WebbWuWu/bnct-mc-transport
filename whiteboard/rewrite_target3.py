# ============================================================
#  rewrite_target3.py —— 默写 sample_target（第 3 次）
#  2026-10-18 Claude 开文件，只写了这段注释；下面全部由你写
# ============================================================
#  要写：
#    1. import（要哪些自己想；函数本身用不到的就不写）
#    2. path_A_prob
#    3. sample_target
#  不写：assert、if __name__ == "__main__": 那一行和下面整块、print
#  不看：lesson07_freegas.py、rewrite_target1.py、rewrite_target2.py、rewrite_sample.py、讲义
#  写完：存盘，把文件名发给 Claude 对拍（同种子跑原函数和你的函数，逐位比）
# ============================================================
import math
from rng import lcg
def path_A_prob(y):
    return y*math.sqrt(math.pi)/(y*math.sqrt(math.pi)+2)
def sample_target(y,rng):
    for k in range(1000):
        xi1=rng.random()
        xi2=rng.random()
        xi3=rng.random()
        if rng.random()<path_A_prob(y):
            x2=-math.log(1-xi1)-math.log(1-xi2)*math.cos((1-xi3)/2*math.pi)**2
        else:
            x2=-math.log((1-xi1)*(1-xi2))
        x=math.sqrt(x2)
        mu_T=2*rng.random()-1
        m=x2+y**2-2*x*y*mu_T
        if rng.random()<math.sqrt(max(m,0.0))/(x+y):
            return x,mu_T,k+1
