# ============================================================
#  rewrite_scatter1.py —— 默写 scatter + alpha（第 1 次，10/23）
#  2026-10-23 Claude 开文件，只写了这段注释；下面全部由你写
# ============================================================
#  要写：
#    1. import（要哪些自己想；函数本身用不到的就不写）
#    2. scatter(E, A, rng)：质心系各向同性弹性散射一次 → 返回 (散射后能量, 实验室系散射角余弦)
#    3. alpha(A)：一次弹性散射后能量的最小比值
#  不写：assert、if __name__ == "__main__": 那一行和下面整块、print
#  不看：lesson04_energy.py、archive\debug_practice_energy.py、lesson07_check1.py
#  写完：存盘，把文件名发给 Claude 对拍（同一组 E、A、种子跑原函数和你的函数，逐位比）
#  写完顺口答一句（物理，不计分）：A = 16 时 alpha 是多少？它和「撞 O 一次平均丢 11%」是什么关系？
# ============================================================
from rng import lcg
import math

def scatter(E,A,rng):
    mu_cm=2*rng.random()-1
    E_new=E*((A**2+2*mu_cm*A+1)/(A+1)**2)
    mu_lab=(A*mu_cm+1)/math.sqrt(A**2+2*A*mu_cm+1)
    return E_new,mu_lab

def alpha(A):
    return (A-1)**2/(A+1)**2
