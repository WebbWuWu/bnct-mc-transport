# -*- coding: utf-8 -*-
# ============================================================
#  lesson01_error.py  --  第1课 任务C：统计误差估计
# ============================================================
#
# 【为什么这一课最重要】
#   蒙特卡洛给出的每一个数，如果不附带误差，就是废数。
#   MCNP / OpenMC / PHITS 的输出里，每个 tally 后面都跟着一个
#   relative error（相对误差 R）。判读规则是行业惯例：
#       R < 0.05   结果可信
#       R < 0.10   一般可信
#       R > 0.10   不可用
#   这一课你要自己把这个 R 算出来。
#
# 【公式】
#   样本均值        x_bar = (1/N) * sum(x_i)
#   样本方差        S^2   = (1/(N-1)) * sum( (x_i - x_bar)^2 )
#                         = (N/(N-1)) * ( sum(x_i^2)/N - x_bar^2 )   <-- 用这个，只需一次遍历
#   均值的标准误    sigma_mean = S / sqrt(N)
#   相对误差        R = sigma_mean / x_bar
#
#   注意 S 和 sigma_mean 的区别：
#     S          描述"单个样本"有多分散 —— 不随 N 变小
#     sigma_mean 描述"均值"有多不准     —— 按 1/sqrt(N) 变小
#   MC 的收敛靠的是后者。

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


def sample_flight(rng, sigma_t):
    xi = rng.random()
    return -math.log(1 - xi) / sigma_t


SIGMA_T = 1.6658


def run(N, seed):
    """跑 N 个历史，返回 (均值, 标准误, 相对误差)"""
    rng = LCG(seed)

    # 【要点】只用两个累加器就够，不需要把百万个样本都存进列表。
    #   这在 GPU 上至关重要 —— 显存装不下上亿个样本，
    #   但两个累加器只占 16 字节。
    s1 = 0.0          # sum of x
    s2 = 0.0          # sum of x^2

    for i in range(N):
        x = sample_flight(rng, SIGMA_T)
        s1 += x
        s2 += x**2                       # TODO: 累加 x 的平方

    mean = s1 / N

    # TODO: 用上面的公式算样本方差 S^2
    #   提示：var = (N/(N-1)) * ( s2/N - mean*mean )
    var = N/(N-1)*(s2/N-mean**2)

    # TODO: 均值的标准误 = sqrt(var) / sqrt(N)
    sigma_mean = math.sqrt(var)/math.sqrt(N)

    # TODO: 相对误差 R = 标准误 / 均值
    R = sigma_mean/mean

    return mean, sigma_mean, R


# ============================================================
#  检验 1：R 是否随 N 按 1/sqrt(N) 下降
# ============================================================
print("=== 检验1：收敛率 ===")
print("      N        均值        标准误        R           1/sqrt(N)    R/(1/sqrt(N))")
for N in [1000, 10000, 100000, 1000000]:
    mean, sm, R = run(N, seed=2026)
    ideal = 1.0 / math.sqrt(N)
    print("%9d  %.6f  %.3e  %.6f   %.6f    %.4f"
          % (N, mean, sm, R, ideal, R / ideal))

# 【思考】最后一列应该接近 1.0。为什么？
#   因为指数分布的标准差恰好等于它的均值（sigma = 1/Sigma_t = mean），
#   所以 R = (sigma/sqrt(N)) / mean = 1/sqrt(N)。
#   换个分布这个比值就不是 1 了，但 1/sqrt(N) 的规律不变。


# ============================================================
#  检验 2：误差棒是不是真的靠谱
#     跑 20 个不同种子，看有多少次"真值落在 mean +- 1*sigma_mean 内"
#     理论上应该约 68%（正态分布 1 sigma 覆盖率）
# ============================================================
print("\n=== 检验2：误差棒可信度（20个种子，N=100000）===")
theory = 1.0 / SIGMA_T
hit1 = 0
hit2 = 0
seeds = [2026, 1, 7, 12345, 99991, 555, 31337, 2718281, 161803, 42,
         101, 202, 303, 404, 505, 606, 707, 808, 909, 1010]
for seed in seeds:
    mean, sm, R = run(100000, seed)
    # 真值离样本均值差了几个标准误
    n_sigma = (mean - theory) / sm
    if abs(n_sigma) < 1.0:
        hit1 += 1
    if abs(n_sigma) < 2.0:
        hit2 += 1
    print("seed %7d   mean=%.6f   偏离 = %+.2f sigma" % (seed, mean, n_sigma))

print("\n落在 1 sigma 内: %d/20 = %.0f%%   （理论 68%%）" % (hit1, hit1 / 20 * 100))
print("落在 2 sigma 内: %d/20 = %.0f%%   （理论 95%%）" % (hit2, hit2 / 20 * 100))
