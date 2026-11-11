# -*- coding: utf-8 -*-
# ============================================================
#  plot_m3.py  --  M3 定稿图：水模体，S(α,β) + 真实弹性角分布都开，本程序 vs OpenMC
#  2026-11-11 Claude 写（README 第一屏要用；用户 11/9 定「Claude 先给画图代码」）。照 plot_water.py（用户 10/27 写的 M2 图）的样子画，
#  多一条虚线：OpenMC 不开 S(α,β) 时的 D_B —— 一眼看出 S(α,β) 把硼剂量往浅处挪了多少。
#
#  怎么跑：Windows、openmc 环境，VS Code 里打开本文件 Ctrl+F5。只读 outputs\ 里的 json，不跑输运，几秒钟。
#  读哪些文件（都是 11/6 跑好的，和 checks\check_M3_seeds.py SA 用的是同一批）：
#     本程序   water_SA_N20000_s3001 … s3010.json（10 × 2 万）+ water_SA_N100000_s1010 / s1011.json（2 × 10 万）= 40 万条
#     OpenMC   openmc_water_sab_s1 … s10.json（开 S(α,β)，10 × 1000 万）
#     虚线     openmc_water_s1 … s10.json（不开 S(α,β)，10 × 1000 万）
#  出什么：outputs\m3_final.png（上格：kerma 深度曲线；下格：本程序 ÷ OpenMC）
#
#  误差怎么算（和 M3 小结 3.2 一致）：
#     本程序：12 个文件的 Σx、Σx² 直接相加，当成一次 40 万条的运行 → mean_se（按历史算，一条历史内的相关已经在里面）
#     OpenMC：10 个种子的均值；SE = 种子之间的散布 / √10（不用它自报的 SE：层间相关、批次相关都自动算进去）
#  ⚠️ 下格每个点的误差棒只管「这一层」。相邻层是一起涨落的（ρ̄ ≈ 0.3），不能数「有几个点在 1σ 外」来下结论 —— 判据看 check_M3_seeds.py。
# ============================================================
# ---- 路径：本文件在 my_code 的子文件夹里，把 my_code 根目录加进搜索路径（和 plot_water.py 同一段）----
import os, sys
_MY_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # my_code 根目录
for _sub in ("whiteboard", "lessons", "openmc_compare", "checks", ""):
    _p = os.path.join(_MY_CODE, _sub) if _sub else _MY_CODE
    if _p not in sys.path:
        sys.path.insert(0, _p)

import json, math
import matplotlib.pyplot as plt
from mcstat import mean_se
from dose import EV2J, G_PER_KG, D, NL, OUT_DIR
from water import rho_w

MY_SMALL = [("water_SA_N20000_s%d.json" % s, 20000) for s in range(3001, 3011)]
MY_BIG = [("water_SA_N100000_s1010.json", 100000), ("water_SA_N100000_s1011.json", 100000)]
OMC_SAB = ["openmc_water_sab_s%d.json" % i for i in range(1, 11)]
OMC_FREE = ["openmc_water_s%d.json" % i for i in range(1, 11)]
NUCS = ["B10", "H1", "O16"]
DOSES = ["D_B", "D_H", "D_O"]
dz = D / NL
COE = EV2J * G_PER_KG / rho_w / dz          # eV·cm（每源中子、每层）→ Gy·cm²（每单位注量），和 plot_water.py 同一个系数

# [判据] 这张图上的数和已经写进 M3 小结的数对不对得上（两条独立的路：小结的数是 11/6 check_M3_seeds.py SA 算的）
#   抓：文件读漏 / 读重、Σx 和 Σx² 配错、核素配错、分子分母反了、开 S 和不开 S 的 OpenMC 文件拿反
REF_TOTAL = {"B10": 1.0009, "H1": 1.0004, "O16": 1.0004}        # 30 层总量 本程序 ÷ OpenMC（M3 小结 4.1 的 S4 列）
REF_SAB = {3: 1.244, 29: 0.477}                                   # OpenMC 开 ÷ 不开 S 的 D_B：第 3 层 = 1.75 cm，第 29 层 = 14.75 cm（M3 小结 4.2）


def load(name):
    with open(os.path.join(OUT_DIR, name), "r") as f:
        return json.load(f)


def my_layers(files, nuc):
    """本程序：把几个文件的 Σx、Σx² 逐层加起来，再按总历史数算均值和标准误。返回 30 层的 (均值, SE)，单位 Gy·cm²。"""
    su = [0.0] * NL
    sq = [0.0] * NL
    n_all = 0
    for name, n in files:
        d = load(name)
        for k in range(NL):
            su[k] += d["kerma_sum"][nuc][k]
            sq[k] += d["kerma_sq"][nuc][k]
        n_all += n
    mean, se = [], []
    for k in range(NL):
        m, _, s = mean_se(su[k], sq[k], n_all)
        mean.append(m * COE)
        se.append(s * COE)
    return mean, se, n_all


def omc_layers(names, tally):
    """OpenMC：几个种子的均值；SE 用种子之间的散布 / √(种子数)。"""
    runs = [load(name)[tally]["mean"] for name in names]
    n = len(runs)
    mean, se = [], []
    for k in range(NL):
        m = sum(r[k] for r in runs) / n
        sd = math.sqrt(sum((r[k] - m) ** 2 for r in runs) / (n - 1))
        mean.append(m * COE)
        se.append(sd / math.sqrt(n) * COE)
    return mean, se


def ratio(m1, s1, m2, s2):
    r, e = [], []
    for k in range(NL):
        r.append(m1[k] / m2[k])
        e.append(r[k] * math.sqrt((s1[k] / m1[k]) ** 2 + (s2[k] / m2[k]) ** 2))
    return r, e


def plot(my, omc, omc_free_B, rat, n_all):
    deep = [(k + 0.5) * dz for k in range(NL)]
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(7, 7), gridspec_kw={"height_ratios": [3, 2]})
    for j, nuc in enumerate(NUCS):
        c = "C%d" % j
        ax1.plot(deep, omc[nuc][0], color=c, label="%s (%s) OpenMC" % (DOSES[j], nuc))
        ax1.errorbar(deep, my[nuc][0], yerr=my[nuc][1], fmt="o", markersize=3, capsize=2, color=c,
                     label="%s (%s) this code" % (DOSES[j], nuc))
        ax2.errorbar(deep, rat[nuc][0], yerr=rat[nuc][1], fmt="o", markersize=3, capsize=2, color=c,
                     label="%s (%s)" % (DOSES[j], nuc))
    ax1.plot(deep, omc_free_B, color="C0", linestyle="--", linewidth=1, label="D_B (B10) OpenMC, free gas (no S(a,b))")
    ax1.set_yscale("log")
    ax2.axhline(1, color="gray", linewidth=0.8)
    ax2.set_ylim(0.8, 1.2)
    ax1.grid(True, which="both", alpha=0.3)
    ax2.grid(True, which="both", alpha=0.3)
    ax1.legend(fontsize=8, ncol=2, loc="lower left")
    ax2.legend(fontsize=8, ncol=3, loc="upper left")
    ax1.set_title("KUR beam, water + 15 ppm B-10, S(a,b) + elastic angular distributions:\n"
                  "this code (%.0e histories) vs OpenMC 0.15.3 (10 seeds x 1e7)" % n_all, fontsize=9)
    ax1.set_ylabel("Kerma per unit fluence (Gy cm^2)")
    ax2.set_ylabel("this code / OpenMC")
    ax2.set_xlabel("Depth (cm)")
    ax2.set_title("neutron kerma only (MT=301), photon dose not included; error bars: 1 SE per layer (layers are correlated)", fontsize=8)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "m3_final.png")
    fig.savefig(path, dpi=150)
    print("图存在", path)


if __name__ == "__main__":
    files = MY_SMALL + MY_BIG
    my, omc, rat = {}, {}, {}
    for nuc in NUCS:
        m1, s1, n_all = my_layers(files, nuc)
        m2, s2 = omc_layers(OMC_SAB, "water %s heating" % nuc)
        my[nuc] = (m1, s1)
        omc[nuc] = (m2, s2)
        rat[nuc] = ratio(m1, s1, m2, s2)
        tot = sum(m1) / sum(m2)
        # 30 层之和 ÷ 30 层之和，和小结里「单独计数的总量」不是同一个计数器，但应是同一个数 → 容差 2e-4（小结只给到 4 位小数）
        assert abs(tot - REF_TOTAL[nuc]) < 2e-4, "%s：30 层总量比 %.5f，M3 小结写的是 %.4f —— 文件读漏 / 读重了？核素配错了？" % (nuc, tot, REF_TOTAL[nuc])
        print("%-3s 30 层总量 本程序 ÷ OpenMC = %.4f（M3 小结 %.4f）" % (nuc, tot, REF_TOTAL[nuc]))
    assert n_all == 400000, "一共应是 40 万条历史，现在是 %d" % n_all
    free_B, _ = omc_layers(OMC_FREE, "water B10 heating")
    for k, want in REF_SAB.items():
        got = omc["B10"][0][k] / free_B[k]
        assert abs(got / want - 1) < 5e-3, "OpenMC 开 ÷ 不开 S 的 D_B 第 %d 层：%.4f，M3 小结写的是 %.3f —— 两组 OpenMC 文件拿反了？" % (k, got, want)
        print("S(α,β) 效应 第 %d 层（%.2f cm）：D_B 开 ÷ 不开 = %.3f（M3 小结 %.3f）" % (k, (k + 0.5) * dz, got, want))
    plot(my, omc, free_B, rat, n_all)
    plt.show()
