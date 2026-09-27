# -*- coding: utf-8 -*-
# ============================================================
#  whiteboard_08_checks.py -- sample_rejection 的五条判据
#  判据由 Claude 写（2026-09-27），被测函数由用户写。
#  运行：python whiteboard_08_checks.py   全过打印「五条全过」
# ============================================================
import math
from rng import lcg
from whiteboard_08_rejection import sample_rejection, Max


class FixedRNG:
    """假随机数发生器：按给定顺序吐数，吐完再要就报错。用于确定性判据。"""
    def __init__(self, seq):
        self.seq = list(seq)
        self.used = 0
    def random(self):
        assert self.used < len(self.seq), "假 rng 的数用完了：函数多要了随机数"
        v = self.seq[self.used]
        self.used += 1
        return v


def f(x):
    return 2.0 * x
def g(x):
    return 1.0
def g_sample(rng):
    return rng.random()
M = 2.0


# ---- 判据 1（确定性）：手算一条完整的抽样过程 ----
# 抓什么：随机数的消耗顺序（先抽 x 再抽 y）、接受条件必须是严格 <、轮数是 k+1 不是 k
# 手算：第 1 轮 x=0.0, y=0.0*2=0.0, f(0)=0 → 0<0 为假 → 拒绝（写成 <= 会在这里错误接受）
#       第 2 轮 x=0.3, y=0.9*2=1.8, f=0.6 → 拒绝
#       第 3 轮 x=0.8, y=0.5*2=1.0, f=1.6 → 接受，应返回 (0.8, 3)
r = FixedRNG([0.0, 0.0, 0.3, 0.9, 0.8, 0.5])
out = sample_rejection(f, g, g_sample, M, r)
assert out == (0.8, 3), "判据1：应返回 (0.8, 3)，实际 %r" % (out,)
assert r.used == 6, "判据1：应恰好用掉 6 个随机数，实际用了 %d 个" % r.used
print("判据1 通过：手算序列返回 (0.8, 3)，用掉 6 个随机数")


# ---- 判据 2（确定性）：撞满上限必须返回 None ----
# 抓什么：for-else 保险丝有没有接上（没接上会返回别的东西，或者死循环）
# 假 rng 永远吐 0.0 → 每轮 x=0, y=0, f=0 → 永远拒绝
r = FixedRNG([0.0] * (2 * Max))
out = sample_rejection(f, g, g_sample, M, r)
assert out is None, "判据2：撞满 %d 轮应返回 None，实际 %r" % (Max, out)
assert r.used == 2 * Max, "判据2：应用掉 %d 个随机数，实际 %d 个" % (2 * Max, r.used)
print("判据2 通过：撞满 %d 轮返回 None" % Max)


# ---- 判据 3：M 给小了必须当场报错 ----
# 抓什么：「f(x) <= M*g(x)」那条 assert 是不是真的接上了、方向对不对
fired = False
try:
    rng = lcg(2026)
    for i in range(1000):
        sample_rejection(f, g, g_sample, 1.5, rng)
except AssertionError as e:
    fired = True
    print("判据3 通过：M=1.5 时报错 ->", e)
assert fired, "判据3：M=1.5 < 2 跑了 1000 次都没报错，M 的检查没接上"


# ---- 判据 4：分布卡方，对照独立的路 x=√ξ（CDF = x²）----
# 抓什么：抽出来的分布不是 f(x)=2x（接受条件写反、y 的范围错、g 与 g_sample 对不上……）
# 期望占比由解析 CDF F(x)=x² 给出，与舍选法的代码路径完全无关
N = 100000
K = 10                                   # 10 个等宽箱
rng = lcg(2026)
counts = [0] * K
rounds = [0] * N
for i in range(N):
    x, n = sample_rejection(f, g, g_sample, M, rng)
    counts[min(int(x * K), K - 1)] += 1
    rounds[i] = n
chi2 = 0.0
for j in range(K):
    a, b = j / K, (j + 1) / K
    E = N * (b * b - a * a)
    chi2 += (counts[j] - E) ** 2 / E
CRIT = 21.666                            # 自由度 K-1=9，显著性 0.01 的临界值
assert chi2 < CRIT, "判据4：chi2=%g 超过临界值 %g（dof=9, 0.01）" % (chi2, CRIT)
print("判据4 通过：chi2=%.3f < %.3f（dof=9）" % (chi2, CRIT))


# ---- 判据 5：平均轮数 = M ----
# 抓什么：接受率错了（M、g、f 三者对不上）。理论：轮数服从几何分布，p=1/M，
#         均值 1/p = M = 2，单样本方差 (1-p)/p² = 2，SE = sqrt(2/N)
mean_n = sum(rounds) / N
se_n = math.sqrt(2.0 / N)
dev = (mean_n - M) / se_n
assert abs(dev) < 4.0, "判据5：平均轮数 %g，偏离理论值 2 达 %.2f 个 SE" % (mean_n, dev)
print("判据5 通过：平均轮数 %.5f，偏离 %.2f 个 SE" % (mean_n, dev))

print("五条全过")
