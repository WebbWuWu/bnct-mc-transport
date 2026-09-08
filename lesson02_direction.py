# -*- coding: utf-8 -*-
# ============================================================
#  lesson02_direction.py  --  第2课 任务A：各向同性方向抽样
# ============================================================
#
# 【问题】
#   把粒子从一维推到三维。方向用一个"单位向量" Omega = (u, v, w) 表示，
#   满足 u^2 + v^2 + w^2 = 1。
#   "各向同性"的意思是：方向均匀地分布在单位球面上，没有任何方向被偏爱。
#
# 【几乎所有人第一次都会掉进去的陷阱】
#   直觉做法：theta 均匀抽 [0, pi]，phi 均匀抽 [0, 2pi]，然后
#       w = cos(theta),  u = sin(theta)*cos(phi),  v = sin(theta)*sin(phi)
#   这是错的。下面 wrong_direction() 就是这个错法，你会亲眼看到它错在哪。
#
# 【为什么错 —— 立体角】
#   球面上的面积元（立体角元）是
#       dOmega = sin(theta) d(theta) d(phi)
#   注意有个 sin(theta) 因子。均匀抽 theta 等于给每个 theta 相同的权重，
#   但赤道附近（theta≈pi/2）的球面环带面积大、两极附近面积小，
#   所以均匀抽 theta 会让粒子在两极堆积。
#
#   正确做法：换元 mu = cos(theta)，则 d(mu) = -sin(theta) d(theta)，于是
#       dOmega = d(mu) d(phi)
#   sin 因子被吸收掉了。所以：
#       mu  均匀抽在 [-1, 1]
#       phi 均匀抽在 [0, 2*pi]
#   这就是"抽 mu 不抽 theta"这条铁律的来源。你会在每一本 MC 教材里见到它。

import math

A = 1664525
C = 1013904223
M = 2**32


class LCG:
    def __init__(self, seed):
        self.state = seed

    def next_uint(self):
        self.state = (A * self.state + C) % M
        return self.state

    def random(self):
        return self.next_uint() / M


# ------------------------------------------------------------
#  错误示范（不用改，跑一下看看它错成什么样）
# ------------------------------------------------------------
def wrong_direction(rng):
    theta = math.pi * rng.random()            # 错在这里：均匀抽 theta
    phi = 2.0 * math.pi * rng.random()
    w = math.cos(theta)
    u = math.sin(theta) * math.cos(phi)
    v = math.sin(theta) * math.sin(phi)
    return u, v, w


# ------------------------------------------------------------
#  正确版本 —— 你来写
# ------------------------------------------------------------
def isotropic_direction(rng):
    """返回一个各向同性的单位方向向量 (u, v, w)"""

    # TODO 1: mu 均匀抽在 [-1, 1]
    #   提示：rng.random() 给的是 [0,1)，怎么线性变换到 [-1,1)？
    mu = rng.random()*2.0-1.0

    # TODO 2: phi 均匀抽在 [0, 2*pi)
    phi = 2.0*math.pi*rng.random()

    # TODO 3: 由 mu 和 phi 组装出 (u, v, w)
    #   w 就是 mu（因为 mu = cos(theta)，而 w = cos(theta)）
    #   sin(theta) = sqrt(1 - mu^2)   （theta 在 [0,pi] 内 sin 恒非负，所以取正根）
    #   u = sin(theta) * cos(phi)
    #   v = sin(theta) * sin(phi)
    sin_theta = math.sqrt(1-mu**2)
    u = sin_theta*math.cos(phi)
    v = sin_theta*math.sin(phi)
    w = mu

    return u, v, w


# ============================================================
#  检验
# ============================================================
def check(name, func, N=200000):
    rng = LCG(seed=2026)
    su = sv = sw = 0.0
    su2 = sv2 = sw2 = 0.0
    max_norm_err = 0.0
    # w 分箱：把 [-1,1] 切成 10 段，统计每段落了多少个
    wbins = [0] * 10

    for i in range(N):
        u, v, w = func(rng)
        su += u; sv += v; sw += w
        su2 += u*u; sv2 += v*v; sw2 += w*w
        # 单位向量检验
        n = math.sqrt(u*u + v*v + w*w)
        if abs(n - 1.0) > max_norm_err:
            max_norm_err = abs(n - 1.0)
        k = int((w + 1.0) / 2.0 * 10)
        if k == 10:
            k = 9
        wbins[k] += 1

    print("\n===== %s =====" % name)
    print("  |Omega| 最大偏离 1 的量 = %.2e   （应该 ~1e-16）" % max_norm_err)
    print("  <u> = %+.5f   <v> = %+.5f   <w> = %+.5f    （都应 ~0）"
          % (su/N, sv/N, sw/N))
    print("  <u^2>=%.5f  <v^2>=%.5f  <w^2>=%.5f   （都应 ~0.3333）"
          % (su2/N, sv2/N, sw2/N))
    print("  w 的10个分箱（每箱应约 %d）:" % (N // 10))
    print("   ", wbins)


check("错误做法：均匀抽 theta", wrong_direction)
check("正确做法：均匀抽 mu", isotropic_direction)

# 【看懂输出】
#   两种做法的 <u>,<v>,<w> 都是 0，<u^2> 也都不容易一眼看出问题。
#   真正暴露错误的是最后那行 w 的分箱：
#     - 正确做法：10 个箱子基本一样多（这叫阿基米德"帽盒定理"：
#       球面在任意等高带上的面积相等，所以 w = cos(theta) 是均匀分布的）
#     - 错误做法：两端的箱子明显多于中间 —— 粒子堆在两极
#
#   这条经验要记住：判断抽样对不对，光看均值往往看不出来，
#   一定要看分布本身。
