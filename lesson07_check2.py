# -*- coding: utf-8 -*-
# ============================================================
#  lesson07_check2.py  --  第 7 课判据 2：平衡谱卡方（v4 的交付判据）
#  2026-09-30 Claude 写。用户的代码在 lesson07_freegas.py，本文件不改它。
#  运行：Ctrl+F5（别用 F5）。约 5–10 分钟。全部通过最后打印「判据 2 全过」。
# ============================================================
#
#  怎么判（讲义 §8 判据 2，9/28 改过的版本）：
#    真库 H-1、只留弹性道（吸收关掉）、大厚度平板（t_cut 内漏不出去）、源 1 eV，
#    每条历史在同一时刻 t_cut 记一个能量（run_slab 第 13 个返回值），
#    N = 10 万条，按 u = E/kT 分 10 个等概率箱，和【密度谱】√u·e^(−u) 做卡方。
#    跑两次：t_cut 和 2·t_cut。两次都要过 —— 说明 t_cut 已经长到「忘掉初始能量」。
#
#  这条抓什么：细致平衡被破坏 ——
#    · 截面没展宽（自由气体配常数截面）
#    · 运动学里温度和截面不同源（kT 用错）
#    · 上散射没实现（靶核静止）
#
#  9/30 容器里在故意改坏的实现上试过（合成 H-1 库：截面按讲义 §4.1 闭式展宽，t_cut = 100 / 200）：
#    截面不展宽            N=2万  → 响，χ² ≈ 2950，平均 E = 1.20 kT
#    字典 kT 用成 600K 那项 N=2万  → 响，χ² ≈ 140
#    靶核静止              N=2万  → 响，全部中子堆在最低一箱
#    🔴 |v_rel| 权重漏掉（舍选第一轮就接受）N=2万 χ² = 5.8 / 10.2；N=10万 χ² = 14.0，平均 E = 1.504 ± 0.004 kT
#       → 【不响】。讲义 §6/§8 说判据 2 能抓这个 —— 错（Claude 的锅，9/30 改讲义）。
#       这一类错由判据 5（lesson07_checks.py：接受率、平均轮数、x 平均）抓，它直接测 sample_target。
#    正确代码              N=10万 → χ² = 14.6，平均 E = 1.499 ± 0.004 kT（不误报）
#
#  🔴 它不能证明什么：
#    · 方向对不对。无限介质里只看能量，「来向被换成随机方向」这类错看不见（9/30 第 1 问）。
#    · 质心公式里 A 放错项：H-1 的 A ≈ 1，放错和放对几乎一样（判据 1 的注释写过）。
#    · 表下限（1e-5 eV）以下那一小撮：10 万条里约 23 条掉下去过，
#      外推写错或者直接杀掉，卡方都看不出来 —— 那部分由 sigma_at_ext 的测试保证。
#    · kT 差 0.13%（0.025335 vs 0.025301）：N = 10 万时只有约 0.5 SE，看不见。
#
#  Maxwell 用哪个 kT：🔴 字典里那个（load_kT 从库里读的 0.025301），
#    不是讲义 §2.1 写的 0.025335（那是 Claude 的错，9/30 已改结论）。
#    这不是「两边同源」：kT 是输入参数；被检验的是输运 + 运动学跑出来的谱。
#
#  参数为什么这么取：
#    dens = 6.4606e22（软组织里 H-1 的数密度，第 6 课用的同一个数）
#    d = 1e6 cm、z_src = d/2：几十次碰撞走不出几厘米，等于无限介质。
#    分 10 个等概率箱：每箱期望 1 万个，远大于 5；等概率箱不用担心尾箱太小。
#    卡方阈值：两次检验，总误报率 1% → 每次 0.005，dof = 10 − 1 = 9 → 临界值 23.589
#      （scipy.stats.chi2.ppf(0.995, 9) 算的；dof 不再减 1，因为 kT 不是从数据里拟合的）。
#
#  列表名：lesson07 的 run_slab 把「t_cut 时刻的能量列表」叫 E_cut（9/30 用户定）；
#          lesson06 的 run_slab 里 E_cut 是「能量阈值」参数。本文件两处都会出现，各有注释。
# ============================================================

import math
from xs_lookup import load_nuclide
from lesson07_freegas import run_slab, load_kT

# ---------------- 由用户按讲义 §8.2 定 ----------------
T_CUT = 150.0         # 9/30 用户定（时间单位 s/√E，1 单位 = 0.723 μs）。合成库实测平均碰撞 52.6 次；2·T_CUT 约 102 次
N     = 100000
SEED  = 2026
# ------------------------------------------------------

CHI2_CRIT = 23.589    # dof = 9，每次检验 α = 0.005
NBIN      = 10


def cdf_density(u):
    """密度谱 √u·e^(−u) 的累积分布：P(E/kT ≤ u) = γ(3/2,u)/Γ(3/2)（讲义 §2.2 的公式）。"""
    return math.erf(math.sqrt(u)) - 2.0 / math.sqrt(math.pi) * math.sqrt(u) * math.exp(-u)


def equal_prob_edges(nbin):
    """二分法找 nbin 个等概率箱的边界（u 值），共 nbin−1 个。"""
    edges = []
    for i in range(1, nbin):
        target = i / nbin
        lo, hi = 0.0, 50.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if cdf_density(mid) < target:
                lo = mid
            else:
                hi = mid
        edges.append(0.5 * (lo + hi))
    return edges


def build_H1_scatter_only():
    """真库 H-1，只留弹性道（MT=2）= 把吸收关掉。字典里的 kT 用用户的 load_kT 读。"""
    E, mts, XS, awr, Q = load_nuclide("H1", "294K")
    i2 = mts.index(2)
    i102 = mts.index(102)
    # 10/9 第 12 条：lesson07 的 run_slab 现在 ① 按 "name" 认核素 ② 从 "Q_value" 取 Q
    #   ③ 函数开头对每个核素找 MT=102 的位置（找不到会 ValueError）。所以字典要带 "name"、"Q_value"，
    #   而且 mts 里必须有 102。做法：102 那一道的截面全填 0 —— 吸收仍然是关掉的，物理不变。
    #   随机数流不变的理由：sigma_at 返回 [σ2, 0.0]，总截面 σ2 + 0.0 == σ2 逐位相同；
    #   order 写成 [[0, 1]]，sample_reaction 扫到第 0 名时 acc/total 已经是 1.0，必在第 0 道返回，
    #   吃的随机数个数和原来 [[0]] 一样（一个 xi）。check2 仍应逐位 12.949 / 7.442。
    nuc = {"name": "H1", "E": E, "mts": [2, 102], "XS": [XS[i2], [0.0] * len(E)],
           "Q_value": [Q[i2], Q[i102]], "awr": awr,
           "dens": 6.4606e22, "kT": load_kT("H1", "294K")}
    return nuc, [[0, 1]]


def chi2_test(E_list, kT, edges):
    counts = [0] * NBIN
    for E in E_list:
        u = E / kT
        b = 0
        while b < NBIN - 1 and u > edges[b]:
            b += 1
        counts[b] += 1
    n = len(E_list)
    expect = n / NBIN
    chi2 = sum((c - expect) ** 2 / expect for c in counts)
    return chi2, counts


def run_one(nuc, order, t_cut, edges):
    d = 1.0e6
    r = run_slab(N, SEED, d, [nuc], 1.0, order, z_src=d / 2, t_cut=t_cut)
    # 10/9 第 12 条：run_slab 返回字典，按键名取（原来是 r[3], r[0], r[1], r[2] / r[9], r[10], r[11] / r[12]）
    coll, n_T, n_R, n_A = r["coll_per_hist"], r["frac_T"], r["frac_R"], r["frac_A"]
    nt_cut, n_sd, n_stuck = r["nt_cut"], r["n_sd"], r["n_stuck"]
    E_cut = r["E_cut"]                     # lesson07：t_cut 时刻的能量列表（不是 lesson06 的能量阈值）

    # 前提：每条历史都活到了 t_cut（没漏、没被吸收、没撞满）。不满足，谱就不是同一时刻的快照。
    assert nt_cut == N and n_T == 0 and n_R == 0 and n_A == 0 and n_stuck == 0 and n_sd == 0, \
        "不是每条历史都在 t_cut 被截断：nt_cut=%d n_T=%g n_R=%g n_A=%g n_sd=%d n_stuck=%d" \
        % (nt_cut, n_T, n_R, n_A, n_sd, n_stuck)
    assert len(E_cut) == N, "能量列表长度 %d ≠ 历史数 %d" % (len(E_cut), N)

    kT = nuc["kT"]
    chi2, counts = chi2_test(E_cut, kT, edges)
    mean_u = sum(E_cut) / N / kT
    var_u = sum((e / kT) ** 2 for e in E_cut) / N - mean_u ** 2
    se_u = math.sqrt(var_u / N)
    print("t_cut = %g：平均碰撞 %.1f 次/历史" % (t_cut, coll))
    print("   各箱计数（期望各 %d）：%s" % (N // NBIN, counts))
    print("   χ² = %.3f（临界 %.3f）   平均 E/kT = %.4f ± %.4f（密度谱 1.5，仅供参考）"
          % (chi2, CHI2_CRIT, mean_u, se_u))
    return chi2, coll


if __name__ == "__main__":
    assert T_CUT is not None and T_CUT > 0, "先按讲义 §8.2 在文件开头填 T_CUT（平均碰撞 ≳ 30 次）"

    edges = equal_prob_edges(NBIN)
    # 自检（Claude 自己的代码）：第一个边界应为 0.292187（scipy gamma(1.5).ppf(0.1) 算的）
    assert abs(edges[0] - 0.2921871870775917) < 1e-9, "等概率箱边界算错：%r" % edges[0]

    nuc, order = build_H1_scatter_only()
    print("H-1 纯散射：awr = %.5f   kT = %.6f eV（库里读的）" % (nuc["awr"], nuc["kT"]))

    results = []
    for tc in (T_CUT, 2 * T_CUT):
        chi2, coll = run_one(nuc, order, tc, edges)
        results.append((tc, chi2, coll))

    # 判据 2
    for tc, chi2, coll in results:
        assert coll >= 30, "t_cut=%g 时平均只碰撞 %.1f 次，不够「忘掉初始能量」（讲义 §8.2 要 ≳ 30）" % (tc, coll)
        assert chi2 < CHI2_CRIT, \
            "判据2：t_cut=%g 的谱和密度谱 Maxwell 不符：χ² = %.3f ≥ %.3f" % (tc, chi2, CHI2_CRIT)
    print("判据 2 全过")
