# -*- coding: utf-8 -*-
# ============================================================
#  lesson08_xs_numbers.py  --  第 8 课核数脚本：讲义里标 ⚠️ 的数，用库核一遍
#  2026-10-09 Claude 写。只读库、只打印，不改任何文件，不 import 你的 run_slab。
#  运行：Ctrl+F5。约 10–30 秒（要打开 5 个 .h5）。最后打印「核数脚本跑完」。
#  跑完把终端输出整段贴给 Claude。
# ============================================================
#
#  它核哪几个数（讲义 = 项目文档\原理讲义_第8课_剂量与隐式俘获.md）
#    第 1 部分  五个核素在 [1e-5 eV, 1e4 eV] 里「开着」的道：MT、名字、是不是冗余、Q 值、0.0253 eV 截面
#               → 回答「D_N 用哪个 MT」「(n,γ) 反应率要加哪些 MT=102」「吸收 = 所有 MT≠2 时实际有哪些道」
#    第 2 部分  讲义 ⚠️ 的数，讲义值 vs 库值：
#               N-14 (n,p) 的 Q（§3 表 0.626 MeV）和 0.0253 eV 截面（§3.2 / §7.2 表 1.83 b）
#               H-1 (n,γ) 0.0253 eV 截面（§7.2 表 0.3326 b）
#               B-10 σ_el/σ_t 在 0.0253 eV（§5.1「约 0.001」）
#    第 3 部分  B-10 两分支 800/801：MT=801 的份额、每次 (n,α) 平均留下的能量，随能量怎么变
#               （§3.1「93.64% 是热能点的比例」、§10 近似 #15）
#
#  记号
#    MT        ENDF 反应类型编号（2 弹性，102 (n,γ)，600 (n,p0)，800/801 (n,α0)/(n,α1)）
#    Q         反应能，库里单位 eV；这里打印时换成 MeV
#    E_thr     这个道在能量网格上第一个存了数的点（threshold_idx 那个点的能量）
#    「开着」  E_thr < 1e4 eV —— 第 8 课的源最高 10 keV、中子只会越撞越慢，
#              所以 E_thr ≥ 1e4 eV 的道在第 8 课里永远是 0，不打印
#    σ(E)      微观截面，barn；在能量网格上线性插值（和你 xs_lookup.sigma_at 同一种插值）
#
#  🔴 为什么自己用 h5py 读，而不调你的 load_nuclide：
#     load_nuclide 会跳过冗余道、也不返回 Q 值和阈值。这里正要看冗余道（MT=103、107）
#     和 Q 值，所以直接读。第 2 条 assert 用冗余道 = 子道之和，当「读对了道」的检验。
# ============================================================

import os
import h5py
import numpy as np
from xs_lookup import LIB

T      = "294K"
E_TH   = 0.0253          # eV，热能点（讲义里那些「热截面」都是在这个能量）
E_MAX  = 1.0e4           # eV，第 8 课源谱上限
NUCS   = ["H1", "C12", "N14", "O16", "B10"]


def read_all(nuc):
    """读一个核素的全部道（含冗余道）。
    返回 E（能量网格，numpy 数组）和一个字典 rx：MT -> {label, redundant, Q, thr, xs(铺满网格)}。
    """
    path = os.path.join(LIB, nuc + ".h5")
    assert os.path.exists(path), "文件不在这里：%s" % path
    rx = {}
    with h5py.File(path, "r") as f:
        root = f[nuc]
        E = root["energy"][T][:]
        for name in root["reactions"].keys():
            r = root["reactions"][name]
            mt = int(r.attrs["mt"])
            ds = r[T]["xs"]
            thr = int(ds.attrs["threshold_idx"])
            xs = ds[:]
            # [Claude assert 1] 抓：截面和能量网格对不上（温度混用 / 读错数据集）—— 对不上时后面插值全错位
            assert len(xs) + thr == len(E), \
                "%s MT=%d：len(xs)=%d + thr=%d ≠ len(E)=%d" % (nuc, mt, len(xs), thr, len(E))
            full = np.zeros(len(E))
            full[thr:] = xs
            label = r.attrs["label"]
            if isinstance(label, bytes):
                label = label.decode()
            rx[mt] = {"label": label,
                      "redundant": int(r.attrs["redundant"]),
                      "Q": float(r.attrs["Q_value"]),
                      "thr": thr,
                      "xs": full}
    return E, rx


def sig(E, rx, mt, e):
    """MT 道在能量 e 处的微观截面（barn），线性插值。"""
    return float(np.interp(e, E, rx[mt]["xs"]))


# ============================================================
#  第 1 部分：每个核素在 [1e-5, 1e4] eV 里开着的道
# ============================================================
print("=" * 78)
print("第 1 部分  [1e-5, 1e4] eV 里开着的道（%s）。冗余 = 1 的是别的道之和，load_nuclide 会跳过它" % T)
print("=" * 78)

data = {}
for nuc in NUCS:
    E, rx = read_all(nuc)
    data[nuc] = (E, rx)
    n_closed = 0
    print("\n%s  （网格 %d 点，共 %d 个道）" % (nuc, len(E), len(rx)))
    print("   %5s  %-14s %4s  %14s  %12s  %14s" % ("MT", "label", "冗余", "Q (MeV)", "E_thr (eV)", "σ(0.0253) b"))
    for mt in sorted(rx):
        r = rx[mt]
        e_thr = E[r["thr"]]
        if e_thr >= E_MAX:
            n_closed += 1
            continue
        print("   %5d  %-14s %4d  %14.6f  %12.4e  %14.6g"
              % (mt, r["label"], r["redundant"], r["Q"] / 1e6, e_thr, sig(E, rx, mt, E_TH)))
    print("   （另有 %d 个道的阈值 ≥ 1e4 eV，第 8 课里永远是 0，没打印）" % n_closed)

# ---- 冗余道 = 子道之和：确认读对了道 ----
# [Claude assert 2] 抓：某个 MT 下面挂的是别的道的数据 / 阈值补零补错位。
#   两条路：库里 NJOY 存好的冗余道  vs  本脚本自己把子道加起来。
#   ⚠️ 抓不到 800 和 801 互换（加法不管顺序）—— 那个看第 3 部分「801 份额」是不是约 0.936。
#
#   容差 REL_ENDF = 1e-6 的来历（10/9 第一版写 1e-9，拍脑袋，在真库上误报 —— Claude 的锅）：
#     ENDF 每个数存在 11 个字符里，形如 3.602260+3，只有 7 位有效数字
#     → 每个数自带最多 5e-7 的相对舍入误差。
#     MT107 自己舍入一次（≤ 5e-7），800 和 801 各舍入一次、加起来（≤ 5e-7），合计 ≤ 1e-6。
#     真库实测在 0.0253 eV 约 2.6e-8，远在界内；读错道会差几个数量级。
REL_ENDF = 1e-6
E_b, rx_b = data["B10"]
worst = 0.0
for e in (E_TH, 1.0, 1.0e3):
    s107 = sig(E_b, rx_b, 107, e)
    s_sum = sig(E_b, rx_b, 800, e) + sig(E_b, rx_b, 801, e)
    rel = abs(s_sum / s107 - 1.0)
    worst = max(worst, rel)
    assert rel < REL_ENDF, \
        "B10 在 %g eV：MT107=%.10g 与 MT800+MT801=%.10g 相对差 %.2e，超过 ENDF 7 位有效数字的舍入界 %.0e" \
        % (e, s107, s_sum, rel, REL_ENDF)
print("\n冗余道检验：B10 MT107 = MT800 + MT801（0.0253 / 1 / 1000 eV 三点，最大相对差 %.2e）✓" % worst)

E_n, rx_n = data["N14"]
np_children = [mt for mt in rx_n if 600 <= mt <= 649]
if 103 in rx_n and len(np_children) > 0:
    s103 = sig(E_n, rx_n, 103, E_TH)
    s_ch = sum(sig(E_n, rx_n, mt, E_TH) for mt in np_children)
    rel = abs(s_ch / s103 - 1.0)
    assert rel < REL_ENDF, \
        "N14 在 0.0253 eV：MT103=%.10g 与 MT600–649 之和=%.10g 相对差 %.2e（子道 %s）" % (s103, s_ch, rel, np_children)
    print("冗余道检验：N14 MT103 = MT%s 之和（0.0253 eV，相对差 %.2e）✓" % (np_children, rel))
else:
    print("⚠️ N14 没有同时存 MT103 和 MT600 系列（MT103 在不在：%s；600 系列：%s）—— 把这一行告诉 Claude"
          % (103 in rx_n, np_children))


# ============================================================
#  第 2 部分：讲义 ⚠️ 的数
# ============================================================
print()
print("=" * 78)
print("第 2 部分  讲义值 vs 库值（比值 = 库 / 讲义）")
print("=" * 78)
print("   %-44s %12s %12s %9s" % ("量", "讲义", "库", "比值"))


def line(name, lec, lib):
    print("   %-44s %12.6g %12.6g %9.4f" % (name, lec, lib, lib / lec))


mt_np = 600 if 600 in rx_n else 103
line("N-14 (n,p) Q，MeV（用 MT=%d）" % mt_np, 0.626, rx_n[mt_np]["Q"] / 1e6)
# （冗余道的 Q_value 在这个库里全是 0 —— 第 1 部分实测 —— 所以 Q 只能从叶子道 600 读）
line("N-14 (n,p) σ(0.0253)，b（MT=%d）" % mt_np, 1.83, sig(E_n, rx_n, mt_np, E_TH))

E_h, rx_h = data["H1"]
line("H-1 (n,γ) σ(0.0253)，b（MT=102）", 0.3326, sig(E_h, rx_h, 102, E_TH))
line("H-1 (n,γ) Q，MeV（讲义 §3.4：2.224 MeV）", 2.224, rx_h[102]["Q"] / 1e6)

s_el = sig(E_b, rx_b, 2, E_TH)
s_t = sum(sig(E_b, rx_b, mt, E_TH) for mt in rx_b if rx_b[mt]["redundant"] == 0)
line("B-10 σ_el/σ_t (0.0253)（讲义：约 0.001）", 0.001, s_el / s_t)
print("      └ B-10 在 0.0253 eV：σ_el = %.6g b，σ_t（叶子道之和）= %.6g b" % (s_el, s_t))
line("B-10 Q(MT800)，MeV", 2.789323, rx_b[800]["Q"] / 1e6)
line("B-10 Q(MT801)，MeV", 2.311472, rx_b[801]["Q"] / 1e6)


# ============================================================
#  第 3 部分：B-10 两分支随能量
# ============================================================
print()
print("=" * 78)
print("第 3 部分  B-10 两分支：MT801 份额、每次 (n,α) 平均留下的能量 Q_avg")
print("           Q_avg = (σ800·Q800 + σ801·Q801) / (σ800 + σ801)；讲义常数 2.341863 MeV")
print("=" * 78)
Q8 = rx_b[800]["Q"] / 1e6
Q9 = rx_b[801]["Q"] / 1e6
Q_CONST = 2.341863
print("   %12s  %12s  %12s  %12s  %14s" % ("E (eV)", "σ800 (b)", "σ801 (b)", "801 份额", "Q_avg (MeV)"))
for e in (1e-5, 0.0253, 0.5, 1.0, 10.0, 100.0, 1e3, 1e4, 1e5, 1e6):
    s8 = sig(E_b, rx_b, 800, e)
    s9 = sig(E_b, rx_b, 801, e)
    q = (s8 * Q8 + s9 * Q9) / (s8 + s9)
    print("   %12.4g  %12.6g  %12.6g  %12.6f  %14.6f" % (e, s8, s9, s9 / (s8 + s9), q))

# 在网格上 [1e-5, 1e4] 里扫一遍：Q_avg 离常数最远多少
mask = E_b <= E_MAX
s8g = rx_b[800]["xs"][mask]
s9g = rx_b[801]["xs"][mask]
qg = (s8g * Q8 + s9g * Q9) / (s8g + s9g)
print("\n   网格上 [1e-5, 1e4] eV 共 %d 点：Q_avg 最小 %.6f，最大 %.6f MeV"
      % (int(mask.sum()), qg.min(), qg.max()))
print("   与常数 2.341863 的最大相对差 = %.3e" % (np.max(np.abs(qg / Q_CONST - 1.0))))

print("\n核数脚本跑完")
