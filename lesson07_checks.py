# -*- coding: utf-8 -*-
# ============================================================
#  lesson07_checks.py  --  第 7 课判据 3、判据 5（只测 sample_target）
#  2026-09-29 Claude 写。用户的代码在 lesson07_freegas.py，本文件不改它。
#  运行：Ctrl+F5。全部通过最后打印「判据 3、5 全过」。
# ============================================================
#
#  判据 3（确定性）—— 抓：选 A 路的概率公式写错 / A、B 两块面积对调
#  判据 5（统计）  —— 抓：① 随机数流错位（每轮消耗数不对）
#                         ② 舍选循环写错（轮数不对、被拒后没重抽、μT 范围错、
#                            A/B 路在 sample_target 里用反）
#  判据 3 只证明 path_A_prob 这一个函数对；它在 sample_target 里有没有用反，
#  由判据 5 的接受率和 x 平均去抓（9/29 在故意改坏的版本上试过，见日志）。
# ============================================================

import math
from rng import lcg
from lesson07_freegas import sample_target, path_A_prob

# ---------------- 判据 3：两个非特殊值，和 40 位十进制手算常数比 ----------------
# 常数来源：Decimal 40 位精度算 y√π/(y√π+2)，√π 取 40 位 —— 和代码那一行是两条路
C_100  = 0.9888421116170458987      # y = 100
C_001  = 0.008784419364871403604    # y = 0.01
for y, c in ((100.0, C_100), (0.01, C_001)):
    p = path_A_prob(y)
    assert abs(p - c) < 1e-12, \
        "判据3：path_A_prob(%g) = %.16f，应为 %.16f，差 %.3g" % (y, p, c, abs(p - c))
print("判据 3 通过：path_A_prob(100) 与 path_A_prob(0.01) 与手算常数差 < 1e-12")


# ---------------- 判据 5 ① 随机数消耗：每轮正好 6 个 ----------------
class CountingRng:
    """包一层发生器，数一数被取走了几个随机数（不改变数的值和顺序）"""
    def __init__(self, rng):
        self.rng = rng
        self.n = 0
    def random(self):
        self.n += 1
        return self.rng.random()

PER_ROUND = 6    # 按用户 9/29 的写法：ξ2、ξ3、ξ4、选路、μT、接受，A/B 两路都是 6 个

crng = CountingRng(lcg(2026))
for i in range(20000):
    before = crng.n
    x, mu_T, k = sample_target(0.985, crng)
    used = crng.n - before
    assert used == PER_ROUND * k, \
        "判据5①：第 %d 次调用返回 %d 轮，却取走了 %d 个随机数（应为 %d）" % (i, k, used, PER_ROUND * k)
print("判据 5① 通过：2 万次调用，每次消耗 = 6 × 轮数")


# ---------------- 判据 5 ② 接受率与 x 平均 对理论 ----------------
def vbar(x, y):
    """对 μT 在 [-1,1] 上平均的相对速度（闭式）"""
    if x < y:
        return y + x * x / (3.0 * y)
    return x + y * y / (3.0 * x)

def simpson(f, a, b, n=4000):
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4 if i % 2 else 2) * f(a + i * h)
    return s * h / 3.0

def theory(y):
    """返回 (平均轮数, x 平均)。接受率 = ∫x²e^(-x²)v̄ / ∫x²e^(-x²)(x+y)"""
    num  = simpson(lambda x: x * x * math.exp(-x * x) * vbar(x, y), 0.0, 8.0)
    num3 = simpson(lambda x: x ** 3 * math.exp(-x * x) * vbar(x, y), 0.0, 8.0)
    den  = y * math.sqrt(math.pi) / 4.0 + 0.5
    return den / num, num3 / num

N = 100000
rng = lcg(2027)
print("\n   y     接受率(实测)  平均轮数 实测/理论 (偏离)     x 平均 实测/理论 (偏离)")
for y in (0.3, 0.985, 3.0):
    k_th, x_th = theory(y)
    sk = sk2 = sx = sx2 = 0.0
    for i in range(N):
        x, mu_T, k = sample_target(y, rng)
        sk += k;  sk2 += k * k
        sx += x;  sx2 += x * x
    k_m = sk / N;  k_se = math.sqrt((sk2 / N - k_m ** 2) / (N - 1))
    x_m = sx / N;  x_se = math.sqrt((sx2 / N - x_m ** 2) / (N - 1))
    dk = (k_m - k_th) / k_se
    dx = (x_m - x_th) / x_se
    print("%6.3f   %.4f        %.4f / %.4f (%+.1f SE)   %.4f / %.4f (%+.1f SE)"
          % (y, 1.0 / k_m, k_m, k_th, dk, x_m, x_th, dx))
    assert abs(dk) < 4.0, "判据5②：y=%g 平均轮数 %.4f 偏离理论 %.4f 达 %.1f SE" % (y, k_m, k_th, dk)
    assert abs(dx) < 4.0, "判据5②：y=%g x 平均 %.4f 偏离理论 %.4f 达 %.1f SE" % (y, x_m, x_th, dx)
print("判据 5② 通过：三个 y 下平均轮数、x 平均都在 4 SE 内")

print("\n判据 3、5 全过")
