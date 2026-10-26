# -*- coding: utf-8 -*-
# ============================================================
#  lesson11_rho_mine.py  --  M2 判据 3（方案丙）事后校准：量【你的程序】自己的层间相关 ρ̄，再重判一次形状检验
#  判据归 Claude，2026-10-26 21:30 写 —— 写在 2 万条 10 个种子的结果出来【之前】，下面的做法和门槛就此定死。
#  被测：lesson11_seeds.py 跑出的 outputs\water_N20000_s2001–s2010.json（你的 run_slab，每个 2 万条）
#  运行：Windows Ctrl+F5（只读 json，不 import openmc，几秒钟）。
# ============================================================
#
#  为什么要做（10/26 晚，用户选乙）
#    10/26 把你的 10 万条和 OpenMC【全部 10 个种子】比：通量形状 Q 有 8/10 出下限 13.12（「太整齐」），
#    H1 heating 的 Q 却在 42.7–49.7、贴着上限 52.34。判据 3 用的 ρ̄ 是 OpenMC 量的（通量用 0.2），
#    而 z 的涨落主要来自你的程序（OpenMC 粒子数是你的 100 倍，SE 只有你的约 1/10）。
#    你用隐式俘获，一条历史带着权重穿很多层 → 层和层一起涨落 → 你的 ρ 可能比 OpenMC（真实吸收）大。
#    ρ̄ 用小了 → Q 被压低：通量出下限、H1 heating 的「过」也可能是被压出来的。
#
#  记号
#    M        独立运行数（10 个种子，同一设置，每个 2 万条）
#    m_k      第 k 层：一次运行的均值 = Σx / N；se_k 用 mcstat.mean_se 算（和 check3 同一个函数）
#    s_k      第 k 层：10 次 m_k 的样本标准差（路 2：实际散多开）
#    se_k     第 k 层：10 次自报 se_k 的均方根（路 1：你的程序说自己散多开）
#    R        sqrt(30 层 (s_k / se_k)² 的平均)，SE 对的话 ≈ 1
#    ρ̄        30 层两两相关系数（435 对）的平均 —— 和 lesson11_check_rho.py 完全同一种算法
#
#  判据（3 条；每条一行「抓什么」）
#    P1 10 个 json 都在，第 1 层通量两两不同
#       抓：种子没传进去（10 次一模一样 → 散布 = 0，下面全是假的）
#    P2 6 个量的 R 都在 [0.44, 1.62]（门槛和 check_rho W3 同一个，来历同 lesson09_check_se.py）
#       抓：你的程序自报的 SE 本身不对（那样谈不上相关）
#    P3 试响：自报 SE 故意 × 0.5，6 个量的 R 都必须出界
#       抓：P2 的门槛太松
#  只量不判：你的 ρ̄ → 重判用的 ρ（规则同 10/24：往小取到 0.1 的整数倍，< 0 按 0）
#
#  重判（事后校准，不覆盖 10/25 种子 1 的结论；χ² 门槛 [13.12, 52.34] 一个字不动）
#    z_k 的相关 ≈ 你的 ρ（OpenMC 那一半的方差只占约 1%），所以形状检验的 ρ̄ 换成上面「用」的那个数。
#    你的 10 万条（s1010）分别对 OpenMC 原版种子 1–10、各向同性种子 1–10，各判 10 次，报：Q 范围、过了几次。
#    同一张表并排给出 10/24 的 OpenMC ρ̄ 判出来的结果，好对照。
# ============================================================
import json
import math
import os

from mcstat import mean_se
from lesson11_check3 import TABLE, NL, R_MAX, CHI2_995, pick, mine_layers, shape_test

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "outputs")
N_RUN = 20000
SEEDS_MINE = list(range(2001, 2011))
R_LO, R_HI = 0.44, 1.62


def sample_sd(x):
    m = sum(x) / len(x)
    return math.sqrt(sum((v - m) ** 2 for v in x) / (len(x) - 1))


def layers_of(d, sk, qk, sub, n):
    su, sq = pick(d, sk, sub), pick(d, qk, sub)
    m, s = [], []
    for k in range(NL):
        mk, _, sek = mean_se(su[k], sq[k], n)
        m.append(mk)
        s.append(sek)
    return m, s


def ratio_R(ms, ses, f=1.0):
    rs = []
    for k in range(NL):
        s_k = sample_sd([m[k] for m in ms])
        se_k = math.sqrt(sum((f * s[k]) ** 2 for s in ses) / len(ses))
        rs.append(s_k / se_k)
    return math.sqrt(sum(r * r for r in rs) / NL)


def mean_corr(ms):
    M = len(ms)
    cols = [[m[k] for m in ms] for k in range(NL)]
    mus = [sum(c) / M for c in cols]
    sds = [sample_sd(c) for c in cols]
    acc, cnt = 0.0, 0
    for i in range(NL):
        for j in range(i + 1, NL):
            cov = sum((cols[i][t] - mus[i]) * (cols[j][t] - mus[j]) for t in range(M)) / (M - 1)
            acc += cov / (sds[i] * sds[j])
            cnt += 1
    return acc / cnt


def rejudge(mine10, omc_files, row, rho):
    name, sk, qk, sub, t30, ttot, mtk, mtq, coe, _ = row
    m1, s1 = mine_layers(mine10, sk, qk, sub)
    Qs = []
    for fn in omc_files:
        with open(os.path.join(OUT_DIR, fn), "r") as f:
            omc = json.load(f)
        m2 = [x * coe for x in omc[t30]["mean"]]
        s2 = [x * coe for x in omc[t30]["se"]]
        Qs.append(shape_test(m1, s1, m2, s2, rho))
    return Qs


if __name__ == "__main__":
    # ---- P1 ----
    runs = []
    for s in SEEDS_MINE:
        p = os.path.join(OUT_DIR, "water_N%d_s%d.json" % (N_RUN, s))
        assert os.path.exists(p), "P1：%s 不在 —— lesson11_seeds.py 还没跑完？" % p
        with open(p, "r") as f:
            runs.append(json.load(f))
    firsts = [r["track_sum"][0] for r in runs]
    assert len(set(firsts)) == len(firsts), "P1：有两次运行第 1 层通量完全相同 —— 种子没传进去？"
    print("P1 过：10 个种子 %d–%d，各 %d 条，结果两两不同\n" % (SEEDS_MINE[0], SEEDS_MINE[-1], N_RUN))

    # ---- P2 + 量 ρ̄ + P3 ----
    print("%-12s %7s %7s %6s %6s" % ("量", "R", "你的ρ̄", "用", "OMCρ̄"))
    use = {}
    for row in TABLE:
        name, sk, qk, sub = row[0], row[1], row[2], row[3]
        ms, ses = [], []
        for r in runs:
            m, s = layers_of(r, sk, qk, sub, N_RUN)
            ms.append(m)
            ses.append(s)
        R = ratio_R(ms, ses)
        rho = mean_corr(ms)
        use[name] = max(0.0, math.floor(rho * 10) / 10)
        print("%-12s %7.3f %7.3f %6.1f %6.1f" % (name, R, rho, use[name], row[9]))
        assert R_LO <= R <= R_HI, "P2：%s 的 R = %.3f 不在 [%.2f, %.2f] —— 你的程序自报的 SE 不可信" % (name, R, R_LO, R_HI)
        assert not (R_LO <= ratio_R(ms, ses, 0.5) <= R_HI), "P3：%s 的 SE × 0.5，P2 居然没响" % name
    print("P2 过：6 个量你自报的 SE 都可信；P3 试响过：SE × 0.5 → 都出界\n")

    # ---- 重判（事后校准）----
    with open(os.path.join(OUT_DIR, "water_N100000_s1010.json"), "r") as f:
        mine10 = json.load(f)
    lo, hi = CHI2_995[NL - 1]
    print("重判形状检验：你的 10 万条 vs OpenMC 各 10 个种子；门槛 [%.2f, %.2f] 不动（R > 0.2 的层照旧剔除）" % (lo, hi))
    print("%-12s %-10s %18s %7s   %18s %7s" % ("量", "OpenMC", "旧 ρ̄：Q 范围", "过", "你的 ρ̄：Q 范围", "过"))
    for row in TABLE:
        name = row[0]
        for tag, pat in (("原版", "openmc_water_s%d.json"), ("各向同性", "openmc_water_iso_s%d.json")):
            files = [pat % s for s in range(1, 11)]
            old = rejudge(mine10, files, row, row[9])
            new = rejudge(mine10, files, row, use[name])
            fmt = lambda L: ("%6.2f – %6.2f" % (min(x["Q"] for x in L), max(x["Q"] for x in L)), "%d/10" % sum(x["ok"] for x in L))
            a, b = fmt(old), fmt(new)
            print("%-12s %-10s %18s %7s   %18s %7s" % (name, tag, a[0], a[1], b[0], b[1]))
    print("\n🔴 这是事后校准：10/25 种子 1 的结论照样报告；M2 小结两套都写，并写明 ρ̄ 是看了 10/26 的种子扫描之后才换的。")
