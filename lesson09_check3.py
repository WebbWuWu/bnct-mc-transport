# -*- coding: utf-8 -*-
# ============================================================
#  lesson09_check3.py  --  第 9 课判据 3：你的程序 vs OpenMC，真问题逐层比
#  判据归 Claude。10/16 14:50 第一版；10/16 17:10 第二版（加 3-5、支持每边两次独立运行）；
#  10/16 17:20 第三版：里程碑按「合并」那一行判（用户 17:18 定 甲），单对只打印、不判。
#  10/16 20:10 第四版：加 3-4（和你的 lesson09_compare.py 存的 compare_kur.json 逐层对 z）。
#  被测：结果 json（都在 outputs\ 下）
#     你的程序：dose_N100000_s1010.json（第 8 课交付那次）、dose_N100000_s1200.json（10/16 换种子独立再跑）
#     OpenMC  ：openmc_kur.json（种子 1）、openmc_kur_s2.json（种子 2），都是 100 批 × 1 万
#  运行：Windows，Ctrl+F5（只用 math / json / os / random，不要 OpenMC）。几秒。
#        每边第一份必须有；第二份没有时会打印「没找到，跳过」，相关的比较不做（不会悄悄当成过）。
#        全过最后打印「判据 3 全过」。没过时先看打印出来的偏差表，再看报错。
# ============================================================
#
#  这一条在比什么
#    同一块软组织（15 ppm 硼）、同一个 KUR 束流，两个程序各算一遍，30 层逐层比 4 个量：
#      通量（径迹长度）、D_B、D_N、(n,γ)。D_H 不在这里（判据 4，周六上午）。
#    这是里程碑 1 本身（讲义 §9 判据 3）。
#
#  记号
#    m、SE     每层的均值、均值的标准误（都按「每源中子」）
#    z_k       第 k 层：(A − B) / sqrt(SE_A² + SE_B²)    —— 讲义 §8；A、B 两次运行随机数独立
#    n1、n2    30 层里 |z| < 1、|z| < 2 的层数
#    χ²        30 层 z² 之和（自由度 30）
#    Q_B、Q_N  硼、氮每次反应放出的能量（eV）；OpenMC 那边是「反应次数」，乘上它才是 eV
#    「合并」   一边的两次独立运行取平均：m = (m₁ + m₂)/2，SE = sqrt(SE₁² + SE₂²)/2
#
#  两条路
#    路 A：你的 run_slab（自由气体、隐式俘获 + 轮盘赌、径迹长度估计反应率）
#    路 B：OpenMC 0.15.3（自由气体、模拟（不开 survival_biasing）、它自己的估计器）
#    共用的只有截面库和问题定义（这正是要求的：同一个问题）。
#    本文件自己从 json 算 z，不 import 你的任何函数 —— 你的 lesson09_compare.py 是另一条
#    「算 z 的路」，3-4 拿两边的 z 互相核对（等 compare 的框架定了再加）。
#
#  一一对应（先确认在比同一个量，纪律 2；讲义 §6）
#    | 量    | 你的 json 键     | OpenMC tally 名      | OpenMC 那边乘 | 单位（每源中子、每层） |
#    | 通量  | track_sum/sq     | "kur flux"           | 1             | cm                   |
#    | D_B   | Bc_sum/sq        | "kur B10 n,a"        | Q_B           | eV                   |
#    | D_N   | Nc_sum/sq        | "kur N14 n,p"        | Q_N           | eV                   |
#    | (n,γ) | ag_sum/sq        | "kur n,gamma"        | 1             | 次                   |
#    D_B 用 Bc（常数 Q）：OpenMC 的 (n,a) 只给次数，乘常数 Q 和 Bc 是同一个量；Br 和 Bc 只差 2.6e-7（第 8 课）。
#    不用换成 Gy：两边乘同一个系数，z 不变（compare 框架第 1 问 甲，10/16 用户定）。
#    ⚠️ D_B、D_N、(n,γ) 在热区截面都 ∝ 1/v，三个反应率是同一个积分乘不同常数 —— 三行 z 几乎一样，
#       算「一个」独立证据，不是三个（10/16 实测）。
#
#  判据（每条旁边写抓什么）
#    3-1 结构 + 量级：OpenMC json 有 7 个名字、每个 30 层、batches = 100、particles = 10000；
#        4 个量 30 层全 > 0；你的 / OpenMC 30 层总和之比在 [0.7, 1.4]；第 1 层 / 第 30 层谁大两边一致
#        抓：tally 名字拼错 / 拿错 tally；核素名或 score 写错 → OpenMC 悄悄给 0（不报错）；
#            忘了乘 Q（比值 ~1e6）、Q_B 和 Q_N 用反（~3.7 倍）、多除 / 少除 dz = 0.5（正好 2 倍）；网格层顺序反了
#    3-2 里程碑 1（框架第一轮第 1 问 丙）：每个量 30 层同时满足
#          15 ≤ n1 ≤ 25；n2 ≥ 26；max|z| ≤ 4；13.79 ≤ χ² ≤ 53.67（χ²₃₀ 的 99% 区间）
#        🔴 判的是「合并」那一行：你的两次平均 vs OpenMC 两次平均（10/16 17:18 用户定 甲）。
#           门槛一个不改，只是把每边两次独立运行的数据都用上。这是看了数据以后定的 —— 日志里如实写：
#           s1010 单独和两次 OpenMC 比都「太好」（通量 1σ 内 28/30、26/30），3-5 证明两边 SE 都对，
#           s1200 单独比全过 → s1010 是碰巧离真值近。
#           4 个单对照样打印（偏差表附表），不参与判定。只有一边有两份时，判单对。
#        抓：两边物理不一致（max|z|、χ² 偏大、n2 偏少）；也抓误差棒估大了（n1 > 25、χ² 偏小）
#    3-3 判据自己会不会响（纪律 4）：用你的均值当「真值」，两边各加和各自 SE 一样大的高斯噪声，
#        造 200 组「两边都对」的假数据，再造两种错的：
#          ① OpenMC × 1.0132（= 10/13 'wo' 归一化那个坑的大小）  ② OpenMC 的 SE × 10（= SE 当成 σ₁、乘了 √批数）
#        要求：「都对」的 ≥ 80% 过 3-2；两种错的 ≥ 95% 被抓到
#        抓：3-2 太松（抓不到 1% 量级的系统偏差）或太紧（两个都对的程序也过不了）
#        ⚠️ 假数据各层独立；真数据相邻层相关（讲义 §8），所以真数据的 n1、χ² 比假数据晃得更厉害。
#    3-4 你的 lesson09_compare.py（另一条「算 z 的路」，你写、没看本文件）存的 outputs\compare_kur.json
#        和本文件「合并」那一行逐个比：4 个量 × 30 层 z，|差| < 1e-9；n1、n2 相等；zmax、χ² |差| < 1e-9
#        抓：compare 里任何一处算错（下标、合并公式、SE 没乘系数、键前缀 / tally 名拿反）；
#            也反过来抓本文件算错 —— 两条路谁错都会响
#        compare_kur.json 不存在时打印「跳过」，结尾不写 3-4 过
#    3-5 同一个程序跑两次（不同种子）互比，用 3-2 同样的门槛（OpenMC 那边连 3 个 heating 也比）
#        抓：这个程序自己的 SE 算错了 —— 两次互比的 χ² 应该 ≈ 30；
#            χ² 明显偏小（n1 > 25）= SE 估大了（比如历史之间不独立、发生器有相关）；χ² 偏大 = SE 估小了
#        10/16 起因：你的 vs OpenMC 通量 1σ 内 28/30（太好了），要分清是你那一次碰巧离真值近，还是 SE 估大了。
#        OpenMC 两次互比：通量 χ² 29.6、热反应 37.7 —— OpenMC 的 SE 没问题（10/16 实测）。
#
#  没过的时候怎么读
#    先看 3-5：同一程序两次互比过了，SE 就可信；再看 3-2：
#      · 整条曲线同号 1–2 SE：可能是一次整体浮动（相邻层相关）→ 看另一对、看「合并」那行
#      · 只在前几层 / 只在最深几层：想想那一段是什么能量的中子在起作用
#      · 某一层 |z| > 4：几乎肯定是真偏差
#    已知差别（偏差表要列，讲义 §10）：#21 非弹当吸收、#22 N-14 (n,α) 等没进打分 —— 预期都 < 0.1%。
# ============================================================

import math
import json
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "outputs")

# (标签, 文件名, 历史数) —— json 里没存 N，按文件名里的 N 写在这里
MINE_RUNS = [
    ("你的 s1010", "dose_N100000_s1010.json", 100000),
    ("你的 s1200", "dose_N100000_s1200.json", 100000),
]
OMC_RUNS = [
    ("OpenMC 种子1", "openmc_kur.json"),
    ("OpenMC 种子2", "openmc_kur_s2.json"),
]
BATCHES = 100                    # 10/15 框架 C 第 ④ 步
PARTICLES = 10000

NL = 30
DZ = 0.5                         # 15 cm / 30 层

Q_B = 2341569.0                  # eV，B-10 (n,α) 常数 Q：lesson07_freegas.py 第 203 行 Q_val（10/9 核数脚本，801 份额 0.937016）
Q_N = 625876.1                   # eV，N-14 MT 600 的 Q（库里读的；lesson08_check4.py 第 63 行）

# (显示名, 你的 json 键前缀, OpenMC tally 名, OpenMC 那边乘的系数)
QUANT = [
    ("flux",    "track", "kur flux",    1.0),
    ("D_B",     "Bc",    "kur B10 n,a", Q_B),
    ("D_N",     "Nc",    "kur N14 n,p", Q_N),
    ("(n,g)",   "ag",    "kur n,gamma", 1.0),
]
HEATING = ["kur H1 heating", "kur C12 heating", "kur O16 heating"]
ALL_OMC_NAMES = [q[2] for q in QUANT] + HEATING

# 里程碑 1 的门槛（框架第一轮第 1 问 丙，10/11 用户定）
N1_LO, N1_HI = 15, 25
N2_LO = 26
ZMAX = 4.0
CHI2_LO, CHI2_HI = 13.79, 53.67   # χ²₃₀ 的 0.5% 和 99.5% 分位

RATIO_LO, RATIO_HI = 0.7, 1.4     # 3-1 量级：只抓「单位 / 系数 / 拿错 tally」这种几倍、几个数量级的错

N_FAKE = 200
FAKE_SEED = 20261016
BIAS = 1.0132                     # 3-3 ①
SE_BLOWUP = 10.0                  # 3-3 ②：√批数 = √100


def mine_layers(r, key, n):
    """你的 json → 每层 (均值, SE)，第 4 课公式：var = (Σx²/N − m²)·N/(N−1)，SE = √(var/N)。"""
    s = r[key + "_sum"]
    sq = r[key + "_sq"]
    m_list, se_list = [], []
    for k in range(NL):
        m = s[k] / n
        var = (sq[k] / n - m * m) * n / (n - 1)
        m_list.append(m)
        se_list.append(math.sqrt(max(var, 0.0) / n))
    return m_list, se_list


def omc_layers(o, name, coef):
    """OpenMC json → 每层 (均值, SE)，都乘 coef（次数 → eV）。"""
    return ([x * coef for x in o[name]["mean"]],
            [x * coef for x in o[name]["se"]])


def merge(a, b):
    """同一个程序两次独立运行合并：均值取平均，SE = sqrt(SE₁² + SE₂²)/2。"""
    (ma, sa), (mb, sb) = a, b
    return ([(x + y) / 2 for x, y in zip(ma, mb)],
            [math.sqrt(x * x + y * y) / 2 for x, y in zip(sa, sb)])


def z_list(ma, sa, mb, sb):
    zs = []
    for k in range(NL):
        den = math.sqrt(sa[k] ** 2 + sb[k] ** 2)
        assert den > 0.0, "第 %d 层两边 SE 都是 0，z 没法算" % (k + 1)
        zs.append((ma[k] - mb[k]) / den)
    return zs


def summary(zs):
    n1 = sum(1 for z in zs if abs(z) < 1.0)
    n2 = sum(1 for z in zs if abs(z) < 2.0)
    zmax = max(abs(z) for z in zs)
    chi2 = sum(z * z for z in zs)
    return n1, n2, zmax, chi2


def failures(n1, n2, zmax, chi2):
    """3-2 的四个门槛，返回没过的那几条（空列表 = 全过）。"""
    bad = []
    if not (N1_LO <= n1 <= N1_HI):
        bad.append("1σ 内 %d 层（应 %d–%d）" % (n1, N1_LO, N1_HI))
    if n2 < N2_LO:
        bad.append("2σ 内 %d 层（应 ≥ %d）" % (n2, N2_LO))
    if zmax > ZMAX:
        bad.append("max|z| = %.2f（应 ≤ %.0f）" % (zmax, ZMAX))
    if not (CHI2_LO <= chi2 <= CHI2_HI):
        bad.append("χ² = %.1f（应在 [%.2f, %.2f]）" % (chi2, CHI2_LO, CHI2_HI))
    return bad


def compare_block(title, pairs, show_z):
    """pairs = [(显示名, (mA, sA), (mB, sB))]。打印逐层 z（可选）和一张偏差表，返回 {显示名: 没过的条目}。"""
    print("  【%s】" % title)
    bad_all = {}
    rows = []
    for disp, a, b in pairs:
        zs = z_list(a[0], a[1], b[0], b[1])
        n1, n2, zmax, chi2 = summary(zs)
        rows.append((disp, n1, n2, zmax, chi2))
        if show_z:
            print("    %s 逐层 z（10 层一行）：" % disp)
            for i in range(0, NL, 10):
                print("      " + " ".join("%+5.2f" % z for z in zs[i:i + 10]))
        bad = failures(n1, n2, zmax, chi2)
        if bad:
            bad_all[disp] = bad
    print("    %-16s %8s %8s %8s %8s %7s  %s" % ("量", "1σ内", "2σ内", "max|z|", "χ²", "χ²/30", "结论"))
    for disp, n1, n2, zmax, chi2 in rows:
        print("    %-16s %5d/30 %5d/30 %8.2f %8.1f %7.2f  %s"
              % (disp, n1, n2, zmax, chi2, chi2 / NL, "过" if disp not in bad_all else "✗ " + "；".join(bad_all[disp])))
    print()
    return bad_all


def load_runs(runs, kind):
    """读存在的 json；每边第一份必须有，后面的没有就打印跳过。返回 [(标签, 数据, 其余)]。"""
    got = []
    for i, run in enumerate(runs):
        path = os.path.join(OUT_DIR, run[1])
        if not os.path.exists(path):
            assert i > 0, "找不到%s的第一份结果 %s" % (kind, path)
            print("  ⚠️ 没找到 %s（%s），跳过和它有关的比较" % (run[1], run[0]))
            continue
        with open(path, "r") as f:
            got.append((run[0], json.load(f), run[2:]))
    return got


if __name__ == "__main__":
    # ---------------- 读 ----------------
    print("=" * 78)
    mines = load_runs(MINE_RUNS, "你的程序")
    omcs = load_runs(OMC_RUNS, "OpenMC")

    # ---------------- 3-1 结构 + 量级 ----------------
    print("3-1 结构 + 量级（先确认两边在比同一个量）")
    for lab, o, _ in omcs:
        assert o.get("batches") == BATCHES and o.get("particles") == PARTICLES, \
            "%s：batches = %r、particles = %r，应为 %d、%d" % (lab, o.get("batches"), o.get("particles"), BATCHES, PARTICLES)
        for name in ALL_OMC_NAMES:
            assert name in o, "%s 里没有 %r —— 名字和 KUR_TALLIES 对不上？现有的键：%r" % (lab, name, list(o.keys()))
            assert len(o[name]["mean"]) == NL and len(o[name]["se"]) == NL, \
                "%s 的 %r 不是 %d 层（%d / %d）" % (lab, name, NL, len(o[name]["mean"]), len(o[name]["se"]))
    for lab, r, _ in mines:
        assert len(r["track_sum"]) == NL and abs(r["dz"] - DZ) < 1e-12, \
            "%s 不是 30 层 × 0.5 cm（%d 层，dz = %r）" % (lab, len(r["track_sum"]), r["dz"])

    # 每个量：每次运行的 (m, SE)
    M = {}   # M[量][标签] = (m, se)
    O = {}
    for disp, key, oname, coef in QUANT:
        M[disp] = {lab: mine_layers(r, key, extra[0]) for lab, r, extra in mines}
        O[disp] = {lab: omc_layers(o, oname, coef) for lab, o, _ in omcs}

    print("  %-6s %-14s %-14s %10s %10s" % ("量", "你的", "OpenMC", "Σ比值", "SE比(O/你)"))
    for disp, key, oname, coef in QUANT:
        for lo, (mb, sb) in O[disp].items():
            zero = [k + 1 for k in range(NL) if mb[k] <= 0.0 or sb[k] <= 0.0]
            assert not zero, "%s / %s：第 %r 层均值或 SE ≤ 0 —— 核素名 / score 写错时 OpenMC 不报错、悄悄给 0" % (lo, disp, zero)
        for lm, (ma, sa) in M[disp].items():
            for lo, (mb, sb) in O[disp].items():
                ratio = sum(ma) / sum(mb)
                se_ratio = (sum(sb) / sum(mb)) / (sum(sa) / sum(ma))
                print("  %-6s %-14s %-14s %10.5f %10.3f" % (disp, lm, lo, ratio, se_ratio))
                assert RATIO_LO <= ratio <= RATIO_HI, \
                    ("%s（%s / %s）：30 层总和之比 = %.4g，不在 [%.1f, %.1f]。多半不是物理差别，是没在比同一个量：\n"
                     "   ~1e6 或 ~1e-6 → 忘了乘 / 多乘了 Q；~3.7 或 ~0.27 → Q_B、Q_N 用反；正好 2 或 0.5 → dz 多除 / 少除了；\n"
                     "   ~1e5 → 你那边没 ÷N；差几十、几百倍 → 拿错了 tally") % (disp, lm, lo, ratio, RATIO_LO, RATIO_HI)
                assert (ma[0] > ma[-1]) == (mb[0] > mb[-1]), \
                    "%s：第 1 层和第 30 层谁大，%s 和 %s 相反 —— 网格层的顺序反了" % (disp, lm, lo)
    print("  3-1 过")

    # ---------------- 3-5 同一个程序跑两次互比（放在 3-2 前面：先确认误差棒可信）----------------
    print("=" * 78)
    print("3-5 同一个程序、不同种子，两次互比（χ² 应 ≈ 30；查的是这个程序自己的 SE 对不对）")
    bad35 = {}
    if len(omcs) == 2:
        pairs = [(d, O[d][omcs[0][0]], O[d][omcs[1][0]]) for d, _, _, _ in QUANT]
        for name in HEATING:
            pairs.append((name.replace("kur ", ""), omc_layers(omcs[0][1], name, 1.0), omc_layers(omcs[1][1], name, 1.0)))
        for d, b in compare_block("%s vs %s" % (omcs[0][0], omcs[1][0]), pairs, False).items():
            bad35["OpenMC " + d] = b
    else:
        print("  ⚠️ OpenMC 只有一份，跳过 OpenMC 互比")
    if len(mines) == 2:
        pairs = [(d, M[d][mines[0][0]], M[d][mines[1][0]]) for d, _, _, _ in QUANT]
        for d, b in compare_block("%s vs %s" % (mines[0][0], mines[1][0]), pairs, True).items():
            bad35["你的 " + d] = b
    else:
        print("  ⚠️ 你的程序只有一份，跳过你的互比 —— 3-2 里「太好了」的原因分不清")

    # ---------------- 3-2 里程碑 1（先全部打印，最后再判）----------------
    print("=" * 78)
    print("3-2 里程碑 1：z = (你的 − OpenMC) / sqrt(SE你的² + SE_OpenMC²)")
    bad_single = {}
    for lm, _, _ in mines:
        for lo, _, _ in omcs:
            pairs = [(d, M[d][lm], O[d][lo]) for d, _, _, _ in QUANT]
            for d, b in compare_block("附表（不判）%s vs %s" % (lm, lo), pairs, False).items():
                bad_single["%s vs %s：%s" % (lm, lo, d)] = b
    if len(mines) == 2 and len(omcs) == 2:
        pairs = [(d, merge(M[d][mines[0][0]], M[d][mines[1][0]]), merge(O[d][omcs[0][0]], O[d][omcs[1][0]]))
                 for d, _, _, _ in QUANT]
        bad32 = {}
        for d, b in compare_block("🔴 里程碑：合并 = 你的两次平均 vs OpenMC 两次平均", pairs, True).items():
            bad32["合并：%s" % d] = b
    else:
        print("  ⚠️ 不是每边都有两份 → 按单对判（第一份 vs 第一份）")
        bad32 = {k: v for k, v in bad_single.items() if k.startswith("%s vs %s" % (mines[0][0], omcs[0][0]))}

    # ---------------- 3-4 和你的 compare 逐层对 ----------------
    print("=" * 78)
    print("3-4 你的 lesson09_compare.py（compare_kur.json）vs 本文件「合并」那一行")
    done34 = False
    cmp_path = os.path.join(OUT_DIR, "compare_kur.json")
    if len(mines) == 2 and len(omcs) == 2 and os.path.exists(cmp_path):
        with open(cmp_path, "r") as f:
            cmp = json.load(f)
        for disp, key, oname, coef in QUANT:
            a = merge(M[disp][mines[0][0]], M[disp][mines[1][0]])
            b = merge(O[disp][omcs[0][0]], O[disp][omcs[1][0]])
            zs = z_list(a[0], a[1], b[0], b[1])
            n1, n2, zmax, chi2 = summary(zs)
            assert oname in cmp, "compare_kur.json 里没有 %r（外层键应是 OpenMC tally 名）；现有：%r" % (oname, list(cmp.keys()))
            c = cmp[oname]
            for k in ["z", "n1", "n2", "zmax", "chi2"]:
                assert k in c, "compare_kur.json[%r] 里没有键 %r（框架定的是 z / n1 / n2 / zmax / chi2）；现有：%r" % (oname, k, list(c.keys()))
            assert len(c["z"]) == NL, "compare_kur.json[%r][\"z\"] 不是 %d 个" % (oname, NL)
            dz_max = max(abs(c["z"][k] - zs[k]) for k in range(NL))
            k_bad = max(range(NL), key=lambda k: abs(c["z"][k] - zs[k]))
            print("  %-14s z 最大差 %.1e（第 %d 层）  n1 %d/%d  n2 %d/%d  zmax 差 %.1e  χ² 差 %.1e"
                  % (oname, dz_max, k_bad + 1, c["n1"], n1, c["n2"], n2, abs(c["zmax"] - zmax), abs(c["chi2"] - chi2)))
            assert dz_max < 1e-9, "%s：第 %d 层 z 不一样：你的 %.6f vs 判据 %.6f" % (oname, k_bad + 1, c["z"][k_bad], zs[k_bad])
            assert c["n1"] == n1 and c["n2"] == n2, "%s：n1 / n2 不一样（你的 %d / %d，判据 %d / %d）" % (oname, c["n1"], c["n2"], n1, n2)
            assert abs(c["zmax"] - zmax) < 1e-9 and abs(c["chi2"] - chi2) < 1e-9, \
                "%s：zmax / χ² 不一样（你的 %.6f / %.6f，判据 %.6f / %.6f）" % (oname, c["zmax"], c["chi2"], zmax, chi2)
        done34 = True
        print("  3-4 过（两条路逐层一致）")
    else:
        print("  ⚠️ 没找到 compare_kur.json（或不是每边两份），跳过 3-4")

    # ---------------- 3-3 判据自己会不会响 ----------------
    print("=" * 78)
    print("3-3 判据自己会不会响：每个量造 %d 组假数据（种子 %d；SE 用 %s、%s 的）" % (N_FAKE, FAKE_SEED, mines[0][0], omcs[0][0]))
    rng = random.Random(FAKE_SEED)
    print("  %-6s %16s %22s %22s" % ("量", "都对：过的比例", "① ×%.4f：抓到的比例" % BIAS, "② SE×%g：抓到的比例" % SE_BLOWUP))
    for disp, key, oname, coef in QUANT:
        ma, sa = M[disp][mines[0][0]]
        mb, sb = O[disp][omcs[0][0]]
        n_ok = n_bias = n_blow = 0
        for t in range(N_FAKE):
            fa = [ma[k] + rng.gauss(0.0, sa[k]) for k in range(NL)]     # 假「你的」
            fb = [ma[k] + rng.gauss(0.0, sb[k]) for k in range(NL)]     # 假「OpenMC」，真值同一个
            if not failures(*summary(z_list(fa, sa, fb, sb))):
                n_ok += 1
            if failures(*summary(z_list(fa, sa, [x * BIAS for x in fb], [x * BIAS for x in sb]))):
                n_bias += 1
            if failures(*summary(z_list(fa, sa, fb, [x * SE_BLOWUP for x in sb]))):
                n_blow += 1
        p_ok, p_bias, p_blow = n_ok / N_FAKE, n_bias / N_FAKE, n_blow / N_FAKE
        print("  %-6s %16.3f %22.3f %22.3f" % (disp, p_ok, p_bias, p_blow))
        assert p_ok >= 0.80, "%s：两边都对的假数据只有 %.0f%% 过得了 3-2 —— 门槛太紧，会误报（告诉 Claude）" % (disp, 100 * p_ok)
        assert p_bias >= 0.95, "%s：OpenMC × %.4f 只抓到 %.0f%% —— 1%% 的系统偏差看不出来（告诉 Claude）" % (disp, BIAS, 100 * p_bias)
        assert p_blow >= 0.95, "%s：SE × %g 只抓到 %.0f%% —— n1 的上限没起作用（告诉 Claude）" % (disp, SE_BLOWUP, 100 * p_blow)
    print("  3-3 过（判据分得出「都对」和这两种错）")

    # ---------------- 判定 ----------------
    print("=" * 78)
    msg = ""
    if bad35:
        msg += "3-5 同一程序两次互比没过（这个程序的 SE 可能算错，先查它，3-2 的结论要等它）：\n" + \
               "\n".join("   %s：%s" % (d, "；".join(b)) for d, b in bad35.items()) + "\n"
    if bad32:
        msg += "3-2 里程碑 1 没过（按文件头「没过的时候怎么读」查）：\n" + \
               "\n".join("   %s：%s" % (d, "；".join(b)) for d, b in bad32.items())
    if msg:
        raise AssertionError("\n" + msg)
    if bad_single:
        print("（附表里没过的单对，只记录、不判：）")
        for d, b in bad_single.items():
            print("   %s：%s" % (d, "；".join(b)))
    print("判据 3 全过（3-1、3-2 合并、3-3、%s3-5）" % ("3-4、" if done34 else "【3-4 没做】、"))
