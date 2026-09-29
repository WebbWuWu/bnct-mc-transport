# -*- coding: utf-8 -*-
# ============================================================
#  lesson07_check1.py  --  第 7 课判据 1：kT → 0 退化到第 6 课
#  2026-09-29 Claude 写。用户的代码在 lesson07_freegas.py，本文件不改它。
#  运行：Ctrl+F5。全部通过最后打印「判据 1 全过」。
# ============================================================
#
#  怎么判：kT = 1e-12 eV（靶核几乎不动），固定 E 和入射方向，
#          抽 10 万次「sample_target → free_gas_kinematics」，
#          ① E'/E 分 20 箱，和 [α, 1] 上的均匀分布做卡方
#          ② 新方向和旧方向的夹角余弦 μ_LAB 取平均，和闭式比   ← 9/29 用户定：留
#
#  ① 抓：四步矢量变换写错（质心速度公式、CM 里速率没保持、没加回质心速度、
#        速度 ↔ 能量换算的单位不一致）
#  ② 抓：E' 对了、但新方向不是由同一个 LAB 速度给出的（拿了 CM 方向、没归一化、
#        又独立抽了一个方向 —— 9/28 补的近似 #11 又回来了）
#
#  9/29 在 Claude 自己写的对照实现上故意改坏 9 处（不是你的代码）：
#    质心除以 A 而不是 1+A → ①  |  CM 里用了 LAB 速率 → ①  |  没加回质心速度 → ①②
#    E' = ½v² 与 E = v² 单位不一致 → ①  |  返回 CM 方向 → ②  |  另抽一个独立方向 → ②
#    方向没归一化 → 单位向量那条  |  沿 z 轴除零 → ±z 那两组直接报 ZeroDivisionError
#    🔴 质心公式里 A 放错项 (A·v_n + V_T)/(1+A) → 只有 ② 响：E'/E 的分布和正确的一模一样
#       （公式对 A 的位置是对称的），只有方向不同。而且 H-1 那两组也不响（A ≈ 1）。
#
#  🔴 它不能证明什么：
#    · 回归判据。kT→0 时上散射本来就不存在，上散射没实现的代码照样通过。
#    · kT→0 时靶核速度约为 0，凡是「乘在靶核速度上」的错都看不见：
#      两个垂直单位向量造错、V_cm 里 A·V_T 写成 V_T（9/29 试过，不响）。
#      这些由函数里的 assert 抓，你写完第 2 步给你。
#
#  参数为什么这么取（纪律：判据参数不取 1、0 这种特殊值）：
#    E = 0.37 eV；A 取库里 H-1 的 awr 0.99917 和 O-16 的 15.8575（都不是整数）；
#    方向取一般方向 (0.36, −0.48, 0.8)，再加两个沿 ±z 轴的 ——
#    后两个是故意挑的「特殊值」：悬着的问题 ② 说的除零就发生在这里，要让它现形。
#
#  卡方阈值：6 组检验，总误报率压在 1% → 每组 0.01/6，dof = 19 → 临界值 42.198
#            （scipy.stats.chi2.ppf(1 - 0.01/6, 19) 算的）。
#            讲义写「p > 0.01」是按一组写的；6 组都用 0.01，正确代码有约 6% 的机会被误判。
# ============================================================

import math
from rng import lcg
from lesson07_freegas import sample_target, free_gas_kinematics

KT    = 1e-12          # eV，靶核几乎静止（不能取 0：y = √(EA/kT) 会除零）
E_IN  = 0.37           # eV
N     = 100000         # 每组散射次数
M     = 20             # 箱数
CHI2_CRIT = 42.198     # dof = 19，单组 p = 0.01/6

NUCLIDES = (("H-1", 0.99917), ("O-16", 15.8575))
DIRS = (("一般方向", (0.36, -0.48, 0.8)),
        ("+z 轴",    (0.0, 0.0, 1.0)),
        ("-z 轴",    (0.0, 0.0, -1.0)))

def alpha(A):
    return ((A - 1.0) / (A + 1.0)) ** 2

def mu_lab_mean(A):
    """静止靶 + CM 各向同性时 μ_LAB 的平均（闭式，推导见日志 9/29）。
    A ≥ 1：2/(3A)；A < 1：1 − A²/3（H-1 的 awr < 1，要用这一支；A = 1 时两支都是 2/3）"""
    if A >= 1.0:
        return 2.0 / (3.0 * A)
    return 1.0 - A * A / 3.0

rng = lcg(2029)
print("kT = %.0e eV   E = %.2f eV   N = %d / 组   箱数 %d   卡方临界 %.3f"
      % (KT, E_IN, N, M, CHI2_CRIT))
print("\n %-5s %-9s %9s   %-7s  %-24s" % ("核素", "入射方向", "α", "① χ²", "② <μ_LAB> 实测 / 闭式 (偏离)"))

for name, A in NUCLIDES:
    al = alpha(A)
    y = math.sqrt(E_IN * A / KT)
    mu_th = mu_lab_mean(A)
    for dname, (u, v, w) in DIRS:
        counts = [0] * M
        s_mu = s_mu2 = 0.0
        for i in range(N):
            x, mu_T, k = sample_target(y, rng)
            E2, u2, v2, w2 = free_gas_kinematics(E_IN, u, v, w, A, KT, x, mu_T, rng)

            # 新方向必须是单位向量 —— 抓：忘了除以 |v_n'|
            n2 = u2 * u2 + v2 * v2 + w2 * w2
            assert abs(n2 - 1.0) < 1e-9, \
                "判据1：%s %s 第 %d 次，新方向长度² = %.12g，不是 1" % (name, dname, i, n2)

            # ① E'/E 落在 [α, 1] 上（kT = 1e-12 带来的越界约 1e-6，给 1e-4 的余量）
            t = (E2 / E_IN - al) / (1.0 - al)
            assert -1e-4 < t < 1.0 + 1e-4, \
                "判据1①：%s %s 第 %d 次，E'/E = %.8f 不在 [α, 1] = [%.6g, 1] 里" % (name, dname, i, E2 / E_IN, al)
            j = min(max(int(t * M), 0), M - 1)
            counts[j] += 1

            # ② μ_LAB = 新方向 · 旧方向
            mu = u * u2 + v * v2 + w * w2
            s_mu += mu
            s_mu2 += mu * mu

        expected = N / M
        chi2 = sum((c - expected) ** 2 / expected for c in counts)
        mu_m = s_mu / N
        mu_se = math.sqrt((s_mu2 / N - mu_m ** 2) / (N - 1))
        d_mu = (mu_m - mu_th) / mu_se
        print(" %-5s %-9s %9.3g   %7.2f  %.5f / %.5f (%+.1f SE)" % (name, dname, al, chi2, mu_m, mu_th, d_mu))
        assert chi2 < CHI2_CRIT, \
            "判据1①：%s %s，E'/E 的 χ² = %.2f ≥ %.3f，不是 [α,1] 上的均匀分布" % (name, dname, chi2, CHI2_CRIT)
        assert abs(d_mu) < 4.0, \
            "判据1②：%s %s，<μ_LAB> = %.5f 偏离闭式 %.5f 达 %.1f SE" % (name, dname, mu_m, mu_th, d_mu)

print("\n判据 1 全过：两个核素 × 三个方向，E'/E 均匀、<μ_LAB> 与闭式相符")
