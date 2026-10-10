# -*- coding: utf-8 -*-
# ============================================================
#  lesson08_check_dose.py  --  第 8 课判据 3c（Gy 换算）+ 判据 5b（KUR 三群源表）
#  判据归 Claude，2026-10-10 写。被测：lesson08_dose.py 里用户写的 to_gy、kur_table
#  运行：Ctrl+F5。几秒。全过最后打印「判据 3c、5b 全过」。
# ============================================================
#
#  为什么不是「第 6、7 条」：每课判据上限 5 条。
#    3c 是判据 3（单位链）的最后一环：3a 验到「每源中子 eV」，3c 验「eV → Gy·cm²」。
#    5b 是判据 5（源抽样）换了新表：10/8 的 S3、S4 期望值按 1/E 算，对 KUR 表不再适用。
#
#  记号
#    「路 A / 路 B」  一条判据两边的两种算法，彼此不共用代码
#    SE              标准误；二项分布的份额 p 抽 n 次：SE = sqrt(p(1−p)/n)
#
#  ------------------------------------------------------------
#  3c  to_gy：eV/源中子 → Gy·cm²
#  ------------------------------------------------------------
#    路 A：用户的 to_gy(s, sq, N, dz)（走 eV → J、g → kg、÷ ρ·dz）
#    路 B：讲义 §2.1 那条路 —— 先 eV → MeV（÷1e6），再 ÷(ρ·dz) 得 MeV·cm²/g，
#          再 × 1.602176634e-10（MeV/g → Gy）。两条路的因子拆法不同，只在最后结果上相遇。
#    测试数据：手编的 4 层 s、sq，N = 12345、dz = 0.437（都不取 1、0.5 这种特殊值，纪律 5）。
#    容差 1e-12（相对）：只有几次乘除，舍入误差 ~1e-16，留足余量。
#    抓：漏 g→kg（差 1000 倍）、×dz 当 ÷dz、1.06 写成换算常数（10/10 你犯过一次）、
#        SE 没乘系数、均值没除 N。
#
#  ------------------------------------------------------------
#  5b  KUR 三群源表
#  ------------------------------------------------------------
#    用 kur_table() 给的表 + 你的 make_cdf / sample_E_table 抽 20 万个能量，数：
#      (1) 三群各占多少  vs  KUR 原文通量比 3.0 : 73 : 4.7（路 B 直接从原文数字算，不走 kur_table）
#      (2) 超热群里 < 10 eV 的占比  vs  1/E 的 ln(10/0.5)/ln(1e4/0.5) = 0.302494
#      (3) 快群里 < 100 keV 的占比  vs  ln(1e5/1e4)/ln(2e6/1e4) = 0.434587
#    每项 |z| < 4。
#    抓：热 / 快份额对调、超热段没乘超热份额、快上限写错（2e7 → (3) 变成 0.302）、
#        热下限写错、箱边界和份额错位一格。
# ============================================================

import math
from rng import lcg
from lesson07_freegas import make_cdf, sample_E_table
from lesson08_dose import to_gy, kur_table


# ============================================================
#  3c
# ============================================================
def check3c():
    N, dz, RHO = 12345, 0.437, 1.06
    s = [8920.7 * N, 1837.0 * N, 321.75 * N, 0.0237 * N]          # 各历史之和（eV）
    sq = [x * x / N * 1.37 for x in s]                              # 平方和：让方差为正、不为 0
    m_a, se_a = to_gy(s, sq, N, dz)
    print("判据 3c  to_gy：eV/源中子 → Gy·cm²（N = %d，dz = %.3f cm，ρ = %.2f）" % (N, dz, RHO))
    worst = 0.0
    for k in range(len(s)):
        mean_ev = s[k] / N
        var_ev = (sq[k] / N - mean_ev ** 2) * N / (N - 1)
        se_ev = math.sqrt(var_ev / N)
        f = 1.0e-6 / (RHO * dz) * 1.602176634e-10          # 路 B：eV → MeV → MeV·cm²/g → Gy·cm²
        m_b, se_b = mean_ev * f, se_ev * f
        r_m = abs(m_a[k] / m_b - 1.0)
        r_s = abs(se_a[k] / se_b - 1.0)
        worst = max(worst, r_m, r_s)
        print("   第 %d 层  均值 %.6e vs %.6e   SE %.6e vs %.6e" % (k, m_a[k], m_b, se_a[k], se_b))
        # [Claude assert] 抓：换算系数里任何一个因子错 —— 都是倍数级，远超 1e-12
        assert r_m < 1e-12 and r_s < 1e-12, \
            "判据3c 第 %d 层：to_gy 给 %.6e ± %.6e，§2.1 那条路给 %.6e ± %.6e（均值差 %.2e、SE 差 %.2e，相对）" \
            % (k, m_a[k], se_a[k], m_b, se_b, r_m, r_s)
    print("   最大相对差 %.1e < 1e-12  ✓" % worst)


# ============================================================
#  5b
# ============================================================
def check5b(n=200000):
    edges, shares = kur_table()
    cdf = make_cdf(shares)
    rng = lcg(20261010)
    n_th = n_epi = n_fast = n_epi_lo = n_fast_lo = 0
    for i in range(n):
        E = sample_E_table(edges, cdf, rng)
        if E < 0.5:
            n_th += 1
        elif E < 1.0e4:
            n_epi += 1
            if E < 10.0:
                n_epi_lo += 1
        else:
            n_fast += 1
            if E < 1.0e5:
                n_fast_lo += 1
    tot = 3.0e7 + 7.3e8 + 4.7e7                       # 路 B：直接用原文通量（Fujiwara 2013 表 1）
    items = [
        ("热群份额 (<0.5 eV)", n_th, n, 3.0e7 / tot),
        ("超热群份额", n_epi, n, 7.3e8 / tot),
        ("快群份额 (>10 keV)", n_fast, n, 4.7e7 / tot),
        ("超热群里 <10 eV", n_epi_lo, n_epi, math.log(10 / 0.5) / math.log(1e4 / 0.5)),
        ("快群里 <100 keV", n_fast_lo, n_fast, math.log(1e5 / 1e4) / math.log(2e6 / 1e4)),
    ]
    print("判据 5b  KUR 三群源表（抽 %d 个）" % n)
    print("   edges 头尾 = %.3g eV … %.3g eV，共 %d 箱" % (edges[0], edges[-1], len(shares)))
    for name, c, m, p in items:
        got = c / m
        se = math.sqrt(p * (1 - p) / m)
        z = (got - p) / se
        print("   %-18s 抽到 %.5f   应为 %.5f   %+.2f SE" % (name, got, p, z))
        # [Claude assert] 抓：份额 / 箱边界 / 段内形状任何一处写错 —— 偏离几十到上千 SE
        assert abs(z) < 4.0, "判据5b %s：抽到 %.5f，应为 %.5f，偏 %.1f SE" % (name, got, p, z)
    # [Claude assert] 抓：热下限 / 快上限没按 10/10 的决定（0.01 eV、2 MeV）
    assert abs(edges[0] - 0.01) < 1e-12 and abs(edges[-1] - 2.0e6) < 1e-6, \
        "表的两端是 %r eV 和 %r eV，10/10 定的是 0.01 eV 和 2e6 eV" % (edges[0], edges[-1])
    print("   五项 |z| < 4，两端 0.01 eV / 2 MeV  ✓")


if __name__ == "__main__":
    check3c()
    print()
    check5b()
    print()
    print("判据 3c、5b 全过")

# ------------------------------------------------------------
#  10/10 容器里在故意改坏的参考实现上试过（每次只改一处）
# ------------------------------------------------------------
#  改坏                                   响在哪   多少
#  to_gy 用 ×dz 代替 ÷dz                  3c       相对差 0.81
#  to_gy 漏 g→kg（×1e3）                  3c       相对差 0.999
#  to_gy 换算常数写成 1.06e-19            3c       相对差 0.34
#  to_gy 的 SE 没乘系数                   3c       SE 差 2.9e15 倍
#  kur_table 热 / 快通量对调              5b       热群 +48.9 SE
#  kur_table 超热段没乘超热份额（再归一） 5b       热群 −8.4 SE
#  kur_table 快上限写成 2e7               5b       快群里 <100 keV −28.7 SE
#  正确实现                               3c 相对差 0；5b 五项都 < 4 SE —— 不误报
