# -*- coding: utf-8 -*-
# ============================================================
#  lesson07_freegas.py  --  第 7 课：自由气体靶核热运动 + 上散射
#  2026-09-28 开文件（Claude：只有这段头注释、import、两个函数的
#  签名行和步骤注释；函数体、run_slab 的改动都由用户写）
# ============================================================
#
#  框架（9/28 用户定的九个决定）见 进度日志 9/28 条目。
#  lesson06_continuous.py 不动 —— 它是假库逐位对拍的基准。
#
#  讲义：项目文档\原理讲义_第7课_热化与自由气体模型.md
#        §4.2 抽靶核  §4.3 四步矢量变换  §8 判据
# ============================================================

import math
from rng import lcg, isotropic_direction

def path_A_prob(y):
    return y*math.sqrt(math.pi)/(y*math.sqrt(math.pi)+2)
def sample_target(y, rng):
    """抽靶核（讲义 §4.2）。
    输入  y   = 中子速度 / 靶核最可几速度 = sqrt(E*A/kT)
          rng = 调用方传进来的发生器（第 4 问：③）
    输出  x, mu_T, 第几轮接受
    """
    # 1. 入口检查 y（Claude 给的 assert 1）
    assert y > 0 and math.isfinite(y), "y=%r 不是正的有限数（E 或 kT 出了问题）" % (y,)
    # 2. 最多 1000 轮：
    for k in range(1000):
        xi2=rng.random()
        xi3=rng.random()
        xi4=rng.random()
        if rng.random()<path_A_prob(y):
            x2=-math.log(1-xi2)-math.log(1-xi3)*math.cos((1-xi4)/2*math.pi)**2
             #A路配方
        else:
            x2=-math.log((1-xi2)*(1-xi3))
        x=math.sqrt(x2)
        mu_T=2*rng.random()-1
        arg = x2 + y**2 - 2*x*y*mu_T
        assert arg > -1e-12, "根号里为负且不是舍入误差：arg=%.3g x=%.6g y=%.6g mu_T=%.6g" % (arg, x, y, mu_T)
        v_rel = math.sqrt(max(arg, 0.0))
        assert v_rel <= x + y + 1e-12, "上界 x+y 失效：x=%.6g y=%.6g mu_T=%.6g v_rel=%.6g" % (x, y, mu_T, v_rel)
        pac=v_rel/(x+y)
        if rng.random()<pac:
            return x,mu_T,k+1
    else:
        assert False, "抽靶核 1000 轮全被拒（正常情况不可能，必是 bug）：y=%r" % (y,)


def free_gas_kinematics(E, u, v, w, A, kT, x, mu_T, rng):
    """四步矢量变换（讲义 §4.3）。
    输出  E', 新方向 u', v', w'
    """
    # 1. 中子速度向量 = 速率 × 方向
    # 2. 靶核速度向量：大小 x × 靶核最可几速度；
    #    与中子方向夹角余弦 = mu_T，绕中子方向的方位角 phi 在 [0, 2pi) 均匀
    #    （❓悬着的 ②：两个与中子方向垂直的单位向量怎么造）
    # 3. 质心速度 = (中子速度 + A × 靶核速度) / (1 + A)
    # 4. CM 里中子速度 = 中子速度 − 质心速度，记下速率
    # 5. CM 里各向同性抽新方向，乘同一个速率
    # 6. 加回质心速度 → LAB 新速度 → E' 和新方向（单位向量）
    pass


# ------------------------------------------------------------
#  Claude 给的三条 assert（成品，去掉行首 # 粘到对应位置）
# ------------------------------------------------------------
# 入口 —— 抓：E 变负或 kT=0 导致 y 是 nan/inf/0
# assert y > 0 and math.isfinite(y), "y=%r 不是正的有限数（E 或 kT 出了问题）" % (y,)
#
# 每轮算完 v_rel 后 —— 抓：mu_T 抽的范围错了（phi 和 mu 那对易混）
# assert v_rel <= x + y + 1e-12, "上界 x+y 失效：x=%.6g y=%.6g mu_T=%.6g v_rel=%.6g" % (x, y, mu_T, v_rel)
#
# for 循环的 else 里 —— 抓：接受概率是 nan 等导致永远不接受（正常撞满概率约 1e-513）
# assert False, "抽靶核 1000 轮全被拒（正常情况不可能，必是 bug）：y=%r" % (y,)


if __name__ == "__main__":
    pass
