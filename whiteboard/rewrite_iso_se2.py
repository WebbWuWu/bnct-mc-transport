# ============================================================
#  rewrite_iso_se2.py —— 默写 lcg / isotropic_direction / mean_se（第 2 次）
#  2026-10-11 Claude 开文件，只写了这段注释；下面全部由你写
# ============================================================
#  要写：
#    1. import（要哪些自己想）
#    2. lcg（发生器，连同它里面的几个方法）
#    3. isotropic_direction
#    4. mean_se
#  不写：assert、if __name__ == "__main__": 那一行和下面整块、print
#  不看：rng.py、mcstat.py、rewrite_iso_se1.py、讲义
#  写完：存盘，把文件名发给 Claude 对拍（同种子比 lcg 前几个数和 isotropic_direction；
#        同一组 s、sq、n 比 mean_se，都要逐位相同）
# ============================================================
import math
a=1664525
c=1013904223
m=2**32
class lcg:
    def __init__(self,seed):
        self.state=seed
    def uint(self):
        self.state=(a*self.state+c)%m
        return self.state
    def random(self):
        return self.uint()/m

def isotropic_direction(rng):
    mu=2*rng.random()-1
    w=mu
    phi=2*rng.random()*math.pi
    sin_theta=math.sqrt(1-mu**2)
    u=sin_theta*math.cos(phi)
    v=sin_theta*math.sin(phi)
    return u,v,w

def mean_se(s,sq,n):
    mean=s/n
    var=n/(n-1)*(sq/n-mean**2)
    se=math.sqrt(var/n)
    return mean,var,se
