# -*- coding: utf-8 -*-
# ============================================================
#  lesson02_slab.py  --  第2课 任务B：你的第一个真正的中子输运程序
# ============================================================
#
# 【问题设置】
#   一块厚度 d 的均匀平板，在 x、y 方向无限大，z 方向从 0 到 d。
#   中子从 z=0 处垂直射入（初始方向 +z，即 (u,v,w) = (0,0,1)）。
#
#   每个中子的命运只有三种，互斥且穷尽：
#       T (透射) : 从 z > d 一侧飞出
#       R (反射) : 被散射回来，从 z < 0 一侧飞出
#       A (吸收) : 在板内被俘获
#   所以必然有  T + R + A = 1 。这是本程序的第一条硬检验。
#
# 【单条历史的流程】—— 就是第0课讲的那个循环
#   1. 出生：z = 0, 方向 = (0,0,1)
#   2. 抽飞行距离 s = -ln(1-xi)/Sigma_t
#   3. 前进：z = z + s*w        （只有 z 方向要跟踪，因为平板在 x,y 无限大）
#   4. 出界判断：z<0 记 R 结束；z>d 记 T 结束
#   5. 没出界就是发生了碰撞，抽反应类型：
#        以概率 Sigma_s/Sigma_t 散射 -> 抽一个新的各向同性方向，回到第 2 步
#        否则                  吸收 -> 记 A，结束
#
# 【为什么 z += s*w】
#   粒子沿单位方向 Omega 走了 s，位移是 s*Omega，它的 z 分量就是 s*w。
#   x、y 分量我们不关心（平板无限大），所以不用跟踪。

import math

A_ = 1664525
C_ = 1013904223
M_ = 2**32


class LCG:
    def __init__(self, seed):
        self.state = seed

    def next_uint(self):
        self.state = (A_ * self.state + C_) % M_
        return self.state

    def random(self):
        return self.next_uint() / M_


def isotropic_direction(rng):
    """把你任务A写好的正确版本抄过来（自己敲，别复制）"""
    mu = rng.random()*2.0-1.0
    phi = 2.0*math.pi*rng.random()
    sin_theta = math.sqrt(1-mu**2)
    u = sin_theta*math.cos(phi)
    v = sin_theta*math.sin(phi)
    w = mu
    return u, v, w


MAX_EVENTS = 10000     # 防止极端情况下死循环（比如 w 恰好接近 0）


def run_slab(N, seed, d, sigma_t, sigma_s):
    """跑 N 条历史，返回 (T, R, A) 三个概率"""
    rng = LCG(seed)

    n_T = 0
    n_R = 0
    n_A = 0
    n_stuck = 0        # 触顶 MAX_EVENTS 的历史数，正常应该是 0

    for i in range(N):
        # ---- 1. 出生 ----
        z = 0.0
        w = 1.0                     # 初始方向 +z，所以 w=1

        for step in range(MAX_EVENTS):
            # ---- 2. 抽飞行距离 ----
            # TODO 1: 用逆变换公式，和第1课一样
            s = -math.log(1-rng.random())/sigma_t

            # ---- 3. 前进 ----
            # TODO 2: 更新 z
            z = z+s*w

            # ---- 4. 出界判断 ----
            # TODO 3: z 小于 0 → 反射出去了
            if z<0:
                n_R += 1
                break
            # TODO 4: z 大于 d → 透射出去了
            if z>d:
                n_T += 1
                break

            # ---- 5. 碰撞：散射还是吸收？----
            # TODO 5: 以 sigma_s/sigma_t 的概率散射
            #   写法：抽一个 [0,1) 随机数，小于这个比例就散射
            if rng.random()<sigma_s/sigma_t:
                # 散射：抽一个新的各向同性方向，只有 w 分量对我们有用
                u, v, w = isotropic_direction(rng)
            else:
                # 吸收：历史结束
                n_A += 1
                break
        else:
            # 【语法点】for...else：循环正常跑完（没被 break）才执行 else
            #   到这里说明碰了 MAX_EVENTS 次还没结束，异常情况
            n_stuck += 1

    if n_stuck > 0:
        print("  !! 警告：有 %d 条历史触顶 MAX_EVENTS" % n_stuck)

    return n_T / N, n_R / N, n_A / N


# ============================================================
#  三个测试
# ============================================================
N = 1000000

print("=== case 1: d=2.0  St=1.0  Ss=0.0  纯吸收（有解析解）===")
T, R, Ab = run_slab(N, 2026, 2.0, 1.0, 0.0)
print("  T = %.6f   R = %.6f   A = %.6f   T+R+A = %.6f" % (T, R, Ab, T+R+Ab))
print("  解析解 exp(-St*d) = %.6f" % math.exp(-1.0*2.0))
print("  说明：没有散射时中子只能直着走，透射率就是 Beer-Lambert 衰减。")
print("        R 必须严格等于 0 —— 不散射就不可能往回走。")

print("\n=== case 2: d=2.0  St=1.0  Ss=0.8  各向同性散射 ===")
T, R, Ab = run_slab(N, 2026, 2.0, 1.0, 0.8)
print("  T = %.6f   R = %.6f   A = %.6f   T+R+A = %.6f" % (T, R, Ab, T+R+Ab))

print("\n=== case 3: d=5.0  St=1.0  Ss=0.9 ===")
T, R, Ab = run_slab(N, 2026, 5.0, 1.0, 0.9)
print("  T = %.6f   R = %.6f   A = %.6f   T+R+A = %.6f" % (T, R, Ab, T+R+Ab))
