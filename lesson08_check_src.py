# -*- coding: utf-8 -*-
# ============================================================
#  lesson08_check_src.py  --  第 8 课判据：源能量表格抽样（1/E 占位表）
#  2026-10-08 Claude 写（判据归 Claude）。被测的三个函数由用户写在 lesson07_freegas.py：
#     table_1overE(edges)          -> 每箱份额（已归一化），长度 = 箱数
#     make_cdf(shares)             -> 累积份额，长度 = 箱数 + 1，首 0、末 1
#     sample_E_table(edges, cdf, rng) -> 一个源能量（eV）
#  运行：Ctrl+F5。全过打印「判据 S1–S4 全过」。
# ============================================================
#
#  记号
#    edges  箱边界，eV。故意取不等宽：等对数宽度时每箱份额全一样，抽错箱也看不出来
#    ln     自然对数
#    SE     标准误 = sqrt(p(1-p)/N)，p 是真值占比，N 是抽样数（不是标准差 σ₁）
#
#  四条判据，两条确定性 + 两条统计。统计的两条各自比的是同一个量：
#    「抽出来的样本里落在某区间的比例」 vs 「1/E 谱在该区间的闭式占比」
#    左边走用户的抽样代码，右边只用 math.log —— 两条独立的路。
# ============================================================

import math
from rng import lcg
from lesson07_freegas import table_1overE, make_cdf, sample_E_table

EDGES = [0.5, 1.0, 3.0, 10.0, 50.0, 200.0, 1000.0, 3000.0, 10000.0]
N     = 100000
SEED  = 20261008
K_SE  = 4.0          # 统计判据允许偏 4 个 SE（正确实现误报概率约 6e-5）

L_TOT = math.log(EDGES[-1] / EDGES[0])     # ln(2e4) = 9.9035

# ---------- S1  份额 ----------
# 抓：份额按 ΔE（每 eV 一样多）而不是按 Δln E（每个对数间隔一样多）算；忘了归一化
shares = table_1overE(EDGES)
assert len(shares) == len(EDGES) - 1, \
    "S1 份额个数 %d ≠ 箱数 %d" % (len(shares), len(EDGES) - 1)
for i in range(len(shares)):
    want = math.log(EDGES[i + 1] / EDGES[i]) / L_TOT
    assert abs(shares[i] - want) < 1e-12, \
        "S1 第 %d 箱 [%g, %g] 份额 %.12f ≠ ln(右/左)/ln(总) = %.12f" \
        % (i, EDGES[i], EDGES[i + 1], shares[i], want)
print("S1 份额 ✓  （第 0 箱 %.6f，第 7 箱 %.6f）" % (shares[0], shares[-1]))

# ---------- S2  CDF ----------
# 抓：CDF 少了开头的 0（和箱边界错位一格 = 栅栏桩）；末尾不是 1（归一化漏了）
cdf = make_cdf(shares)
assert len(cdf) == len(EDGES), \
    "S2 CDF 长度 %d ≠ 箱边界个数 %d（应 = 箱数 + 1，开头那个 0 在不在？）" % (len(cdf), len(EDGES))
assert cdf[0] == 0.0, "S2 CDF 第一个数应为 0，实际 %r" % (cdf[0],)
assert abs(cdf[-1] - 1.0) < 1e-12, "S2 CDF 最后一个数应为 1，实际 %.15f" % cdf[-1]
for i in range(len(shares)):
    assert abs((cdf[i + 1] - cdf[i]) - shares[i]) < 1e-12, \
        "S2 cdf[%d]-cdf[%d] = %.12f ≠ 第 %d 箱份额 %.12f" % (i + 1, i, cdf[i + 1] - cdf[i], i, shares[i])
print("S2 CDF  ✓")

# ---------- 抽 N 个 ----------
rng = lcg(SEED)
n_lo = 0          # 落在 [0.5, 10) eV
n_2k = 0          # 落在 [0.5, 2000) eV
for _ in range(N):
    E = sample_E_table(EDGES, cdf, rng)
    assert EDGES[0] <= E <= EDGES[-1], "抽出的能量 %r 不在 [%g, %g] 里" % (E, EDGES[0], EDGES[-1])
    if E < 10.0:
        n_lo += 1
    if E < 2000.0:
        n_2k += 1

# ---------- S3  0.5–10 eV 占比 ----------
# 抓：抽箱抽错（CDF 错位一格、比较号写反、用错箱号去取边界）。
# 10 eV 正好是箱边界 —— 所以这条只查「抽哪个箱」，箱内怎么抽它看不见（S4 管）。
p3 = math.log(10.0 / 0.5) / L_TOT            # 0.302494
se3 = math.sqrt(p3 * (1 - p3) / N)
f3 = n_lo / N
assert abs(f3 - p3) < K_SE * se3, \
    "S3 0.5–10 eV 占比 %.5f，1/E 真值 %.5f，差 %.1f 个 SE（容许 %g）" % (f3, p3, (f3 - p3) / se3, K_SE)
print("S3 0.5–10 eV  抽样 %.5f  真值 %.5f  差 %+.2f SE ✓" % (f3, p3, (f3 - p3) / se3))

# ---------- S4  < 2 keV 占比 ----------
# 抓：箱内在 E 上均匀抽（应在 ln E 上均匀）；抽箱和箱内用了同一个随机数。
# 2 keV 在箱 [1000, 3000] 的中间：E 上均匀会给 0.8230，ln E 上均匀给 0.8375，差约 12 个 SE。
p4 = math.log(2000.0 / 0.5) / L_TOT          # 0.837494
se4 = math.sqrt(p4 * (1 - p4) / N)
f4 = n_2k / N
assert abs(f4 - p4) < K_SE * se4, \
    "S4 <2 keV 占比 %.5f，1/E 真值 %.5f，差 %.1f 个 SE（容许 %g）" % (f4, p4, (f4 - p4) / se4, K_SE)
print("S4 <2 keV     抽样 %.5f  真值 %.5f  差 %+.2f SE ✓" % (f4, p4, (f4 - p4) / se4))

print("判据 S1–S4 全过")
