# -*- coding: utf-8 -*-
# ============================================================
#  lesson08_checks.py  --  第 8 课判据 1–3（判据归 Claude，2026-10-09 写）
#  被测：lesson07_freegas.py 的 run_slab（用户写；10/9 改动清单 1–12 全部完成后的版本）
#  运行：Ctrl+F5（别用 F5）。判据 3 几秒；判据 2 要跑真库组织两遍（各 2 万条），约 3–10 分钟。
#        全过最后打印「第 8 课判据 1–3 全过」。
# ============================================================
#
#  记号
#    SE      标准误（均值的误差棒）= sqrt(单条历史方差 / N)，不是标准差 σ₁
#    z       偏离几个 SE = (A − B) / sqrt(SE_A² + SE_B²)，A、B 来自两次独立运行
#    Σ       宏观截面 cm⁻¹ = 数密度 × 微观截面(b) × 1e-24
#    「每源中子」 所有返回的 _sum 除以 N
#
#  判据 1  回归：退回第 7 课 —— 就是 lesson07_check2.py。10/9 20:42 用户跑过：χ² = 12.949 / 7.442，
#          平均 E/kT = 1.5043，与 9/30 逐位相同。本文件不重复跑（它要 5–10 分钟），只在最后打印提醒。
#          它抓：改的过程中把第 7 课的东西碰坏了。它【不】证明第 8 课新加的东西对（纪律 4）。
#
#  判据 3  单位链 + Q 值（先跑，几秒）
#    3a  假库、纯吸收（散射截面只留 1e-4 的零头）、单能、垂直入射：第 k 层的
#        硼逐道 / 硼常数 / 氮 / (n,γ) 每源中子的值都有闭式：
#            系数 / Σ_t × ( e^(−Σ_t·z_k) − e^(−Σ_t·z_(k+1)) )
#        系数 = 数密度 × 1e-24 × Σ(σ_道 × Q_道)（(n,γ) 不乘 Q）。
#        两条路：run_slab 一步一步记出来的账  vs  纸上的指数公式。
#        抓：BARN 漏乘 / 多乘、Q 乘错道（800↔801）、氮少乘 Q、/N 或层号错位、wt 没乘或乘了两次。
#        参数都不取 1、0（纪律 5）：Σ_t ≈ 0.4137 cm⁻¹、dz = 0.45 cm。
#    3b  Q 值：load_nuclide 读进字典的 Q_value  vs  直接用 h5py 读库的 Q_value 属性（两条路）。
#        抓：Q 列表和 mts 错位（append 位置放错、冗余道没跳过）。
#
#  判据 2  隐式 vs 真吸收（慢，最后跑）
#    真库软组织（ICRU-44：H 10.2%、C 14.3%、N 3.4%、O 70.8%，ρ = 1.06）+ ¹⁰B 15 ppm 进输运，
#    1/E 源 [0.5 eV, 10 keV]（第 8 课的占位源），表面垂直入射，d = 15 cm，NL = 30。
#    同一问题跑两遍、两个不同种子：implicit=False（真吸收）和 implicit=True（隐式 + 轮盘赌）。
#    五个量逐层比：通量、硼逐道、氮、(n,γ)、D_H。每层 |z| < 4。
#    两条路：两种输运方式（期望相同，随机数流完全不同）。
#    抓：打折因子写错（Σ_a/Σ_t 当成 Σ_s/Σ_t、或打折两次）、轮盘赌赢了提权写错（提到 wt_lo 而不是 wt_hi）、
#        通量 / 剂量没乘 wt、D_H 用了打折前的 wt。
#    🔴 不抓：两种模式【共用】的部分（径迹长度系数、截面查表、运动学）—— 那些错两边一起错，
#        z 照样小。那一半归判据 3。
#    阈值 4 SE 的来历：5 个量 × 30 层 = 150 次比较，单次 |z| ≥ 4 的概率 6.3e-5，
#        全部正确时误报概率约 150 × 6.3e-5 ≈ 1%。
#
#  10/9 容器里在故意改坏的实现上试过（假库组织：五核素、道和截面量级同真库；N = 2 万）：
#    改坏                                  响不响
#    B1 隐式打折两次 (Σ_s/Σ_t)²             判据 2 响：通量第 0 层 −6.9 SE（N=5000 就响）
#    B2 轮盘赌赢了提到 wt_lo 而不是 wt_hi   判据 2 响：通量第 16 层 −4.5 SE（N=5000 时只有 3.3 SE、不响 → 所以 N 取 2 万）
#    B3 D_H 用打折前的 wt                   🔴 不响（N=2 万，最大 2.1 SE）。原因：D_H 几乎全来自 keV 区的 H 散射，
#                                           那里 Σ_s/Σ_t ≈ 0.9999，打折前后差不到 0.1%，这个问题里本来就看不出来。
#                                           「wt 用打折后」这件事靠讲义自测第 5 题的推导保证，不靠这条判据。
#    B4 拆层里通量没乘 wt                   判据 2 响：通量第 0 层 +10.4 SE
#    B5 800 和 801 的 Q 对调                判据 3a 响：硼逐道 +216 SE
#    B6 氮系数少乘 Q                        判据 3a 响（差 6 个数量级）
#    正确代码（今晚用户的版本，N=2 万）     3a 最大 1.05 SE；判据 2 最大 2.09 SE —— 不误报。容器里判据 2 用时 188 s
# ============================================================

import math
import os
import numpy as np
from xs_lookup import load_nuclide, LIB
from lesson07_freegas import run_slab, load_kT, table_1overE

# ------------------------------------------------------------
#  公用：从 _sum / _sq 算每层均值和 SE
# ------------------------------------------------------------
def layer_mean_se(s, sq, N):
    """s、sq 是长 NL 的列表（各历史之和、各历史平方之和）。返回 (均值列表, SE 列表)。"""
    m = [x / N for x in s]
    se = []
    for k in range(len(s)):
        var = (sq[k] / N - m[k] ** 2) * N / (N - 1)
        se.append(math.sqrt(max(var, 0.0) / N))
    return m, se


# ============================================================
#  判据 3a：单位链闭式（假库）
# ============================================================
def check3a():
    E = [1.0e-5, 2.0e7]
    S_TINY = 1.0e-4                       # 散射截面只留零头（run_slab 里 sigma_s 做分母，不能是 0）
    # 假 B-10：800、801 两道常数截面；Q 用库里的真值
    sB800, sB801 = 37.3, 554.9            # b
    nB = 5.27e20                          # cm⁻³
    # 假 N-14：(n,p) 一道
    sN600 = 1.83
    nN = 6.1e21
    # (n,γ)：两个核素各有一点
    sB102, sN102 = 0.211, 0.075
    Q800, Q801, Q600 = 2789323.0, 2311472.0, 625876.0
    B = {"name": "B10", "E": E, "mts": [2, 102, 800, 801],
         "XS": [[S_TINY] * 2, [sB102] * 2, [sB800] * 2, [sB801] * 2],
         "Q_value": [0.0, 11456000.0, Q800, Q801], "awr": None, "dens": nB, "kT": None}
    Nn = {"name": "N14", "E": E, "mts": [2, 102, 600],
          "XS": [[S_TINY] * 2, [sN102] * 2, [sN600] * 2],
          "Q_value": [0.0, 10833390.0, Q600], "awr": None, "dens": nN, "kT": None}
    order = [[3, 2, 1, 0], [2, 1, 0]]

    BARN = 1.0e-24
    Sig_t = nB * BARN * (S_TINY + sB102 + sB800 + sB801) + nN * BARN * (S_TINY + sN102 + sN600)
    c_Br = nB * BARN * (sB800 * Q800 + sB801 * Q801)
    c_Bc = nB * BARN * (sB800 + sB801) * 2341569.0
    c_N = nN * BARN * sN600 * Q600
    c_g = nB * BARN * sB102 + nN * BARN * sN102

    N, NL, d = 100000, 12, 5.4            # dz = 0.45 cm
    r = run_slab(N, 20261009, d, [B, Nn], 0.0253, order, z_src=0.0, implicit=False, NL=NL)
    dz = r["dz"]
    assert abs(dz - 0.45) < 1e-12, "dz=%r，应为 0.45（NL 参数没传进去？）" % (dz,)

    print("判据 3a  单位链闭式（假库：Σ_t = %.4f /cm，dz = %.2f cm，N = %d）" % (Sig_t, dz, N))
    print("   %-8s %4s %14s %14s %8s" % ("量", "层", "程序", "闭式", "偏 SE"))
    worst = 0.0
    for key, c, name in (("Br", c_Br, "硼逐道"), ("Bc", c_Bc, "硼常数"), ("Nc", c_N, "氮"), ("ag", c_g, "(n,γ)")):
        m, se = layer_mean_se(r[key + "_sum"], r[key + "_sq"], N)
        for k in range(NL):
            want = c / Sig_t * (math.exp(-Sig_t * k * dz) - math.exp(-Sig_t * (k + 1) * dz))
            z = (m[k] - want) / se[k]
            worst = max(worst, abs(z))
            if k in (0, 5, 11):
                print("   %-8s %4d %14.6e %14.6e %+8.2f" % (name, k, m[k], want, z))
            # [Claude assert] 抓：单位链上任何一个因子错（BARN、Q、/N、层号、wt）—— 都是倍数级或错层，远超 4 SE
            assert abs(z) < 4.0, \
                "判据3a %s 第 %d 层：程序 %.6e ≠ 闭式 %.6e（偏 %.1f SE）" % (name, k, m[k], want, z)
    print("   全部 4 个量 × %d 层，最大偏离 %.2f SE  ✓" % (NL, worst))


# ============================================================
#  判据 3b：Q 值两条路
# ============================================================
def check3b():
    import h5py
    for nuc, mt_list in (("B10", (800, 801)), ("N14", (600,))):
        E, mts, XS, awr, Q = load_nuclide(nuc, "294K")
        with h5py.File(os.path.join(LIB, nuc + ".h5"), "r") as f:
            for mt in mt_list:
                q_direct = float(f[nuc]["reactions"]["reaction_%03d" % mt].attrs["Q_value"])
                q_dict = Q[mts.index(mt)]
                # [Claude assert] 抓：Q 列表和 mts 错位（第 i 个 Q 不是第 i 个 MT 的）
                assert q_dict == q_direct, \
                    "判据3b %s MT=%d：load_nuclide 给的 Q=%r，库里直接读的 Q=%r" % (nuc, mt, q_dict, q_direct)
                print("判据 3b  %s MT=%d  Q = %.1f eV（两条路逐位相同）✓" % (nuc, mt, q_dict))


# ============================================================
#  判据 2：隐式 vs 真吸收（真库组织）
# ============================================================
M_N = 1.00866491588          # 中子质量，g/mol（awr × M_N = 摩尔质量）
N_A = 6.02214076e23
RHO = 1.06
W = {"H1": 0.102, "C12": 0.143, "N14": 0.034, "O16": 0.708, "B10": 15e-6}
E_KEEP = 2.0e4               # 只留阈值低于 20 keV 的道：源最高 10 keV，更高阈值的道在这一课里恒为 0（省时间）


def build_tissue():
    nucs, orders = [], []
    for nuc in ("H1", "C12", "N14", "O16", "B10"):
        E, mts, XS, awr, Q = load_nuclide(nuc, "294K")
        keep = [i for i in range(len(mts)) if np.any(XS[i][E < E_KEEP] > 0.0)]
        mts_k = [mts[i] for i in keep]
        XS_k = [XS[i] for i in keep]
        Q_k = [Q[i] for i in keep]
        dens = W[nuc] * RHO * N_A / (awr * M_N)
        d = {"name": nuc, "E": E, "mts": mts_k, "XS": XS_k, "Q_value": Q_k,
             "awr": awr, "dens": dens, "kT": load_kT(nuc, "294K")}
        i0 = int(np.searchsorted(E, 0.0253)) - 1
        orders.append(sorted(range(len(mts_k)), key=lambda t: XS_k[t][i0], reverse=True))
        nucs.append(d)
        print("   %-4s 留 %2d 个道（共 %2d），dens = %.4e /cm³" % (nuc, len(mts_k), len(mts), dens))
    return nucs, orders


def check2(N=20000):
    print("判据 2  隐式 vs 真吸收（真库组织 + 15 ppm ¹⁰B，1/E 源，d = 15 cm，NL = 30，各 N = %d）" % N)
    nucs, orders = build_tissue()
    edges = [0.5, 1.0, 3.0, 10.0, 50.0, 200.0, 1000.0, 3000.0, 10000.0]
    src = (edges, table_1overE(edges))
    ra = run_slab(N, 1009, 15.0, nucs, None, orders, z_src=0.0, implicit=False, NL=30, src=src)
    print("   真吸收跑完：平均碰撞 %.1f 次/历史，透射 %.4f，反射 %.4f，吸收 %.4f"
          % (ra["coll_per_hist"], ra["frac_T"], ra["frac_R"], ra["frac_A"]))
    ri = run_slab(N, 2010, 15.0, nucs, None, orders, z_src=0.0, implicit=True, NL=30, src=src)
    print("   隐式跑完：  平均碰撞 %.1f 次/历史，轮盘赌杀掉 %d 条，透射权重 %.4f，反射权重 %.4f"
          % (ri["coll_per_hist"], ri["n_ki"], ri["nt_w"] / N, ri["nr_w"] / N))

    worst_all = 0.0
    for key, name in (("track", "通量"), ("Br", "硼逐道"), ("Nc", "氮"), ("ag", "(n,γ)"), ("DH", "D_H")):
        ma, sa = layer_mean_se(ra[key + "_sum"], ra[key + "_sq"], N)
        mi, si = layer_mean_se(ri[key + "_sum"], ri[key + "_sq"], N)
        worst, kw = 0.0, -1
        for k in range(30):
            s = math.sqrt(sa[k] ** 2 + si[k] ** 2)
            if s == 0.0:
                continue
            z = (mi[k] - ma[k]) / s
            if abs(z) > worst:
                worst, kw = abs(z), k
            # [Claude assert] 抓：隐式那条路的权重处理错（打折、轮盘赌提权、乘 wt）—— 期望偏离，N 够大必响
            assert abs(z) < 4.0, \
                "判据2 %s 第 %d 层：隐式 %.5e ± %.2e vs 真吸收 %.5e ± %.2e，偏 %.1f SE" \
                % (name, k, mi[k], si[k], ma[k], sa[k], z)
        print("   %-6s 第 1 层：真吸收 %.4e ± %.1e  隐式 %.4e ± %.1e   30 层最大 |z| = %.2f（第 %d 层）"
              % (name, ma[0], sa[0], mi[0], si[0], worst, kw))
        worst_all = max(worst_all, worst)
    print("   五个量 × 30 层，最大 |z| = %.2f < 4  ✓" % worst_all)


if __name__ == "__main__":
    check3a()
    print()
    check3b()
    print()
    check2()
    print()
    print("判据 1（回归）：见 lesson07_check2.py，10/9 20:42 已跑 χ² = 12.949 / 7.442，与 9/30 逐位相同。")
    print("第 8 课判据 1–3 全过")
