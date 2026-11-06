# -*- coding: utf-8 -*-
# ============================================================
#  check3_m3_recal.py  --  M3 判据 3 丙：开 S(α,β) 时用 OpenMC 自比重标形状门槛（11/6 用户定甲）
#  判据归 Claude，2026-11-06 17:35 写 —— 在 OpenMC 种子 2–10 跑完之前写下（门槛规则先定，再看数）。
#  运行：Windows Ctrl+F5（只读 json，几秒钟）。要先在 WSL 跑完 explore/m3_water_openmc.py 2 3 4 5 6 7 8 9 10。
#
#  为什么要重标（照实写，看了数据之后才发现的）：
#    原门槛 = χ²(层数 − 1) 的 [0.5%, 99.5%]，假设 30 层之间「均匀相关」（任意两层相关都是 ρ̄）。
#    开 S(α,β) 后，你的程序两个种子自己比（SA1010 vs SA1011）通量 Q = 70.3、B10 Q = 75.9 —— 同一个程序也判不过；
#    而且 Q 两头都出（13.7、15.6 和 57、76）。→ 原门槛的「零分布」在开 S 时不对：热中子一条历史撞约 150 次、
#    在几厘米里来回扩散，相邻几层一起高一起低、远的层不相关，不是均匀相关。
#  怎么重标（同 10/23 M1 重标 D_H 门槛的做法）：
#    OpenMC 开 S 的 10 个种子（各 1000 万）两两比，45 对，每对用同一个公式、同一个 ρ̄ 算 Q（物理一样 → 层间相关结构一样）。
#    新门槛 = [45 对里的最小 Q, 45 对里的最大 Q]（10/23 规则：上限取自比 45 对的最大值）。只改 C3 的门槛，C1、C2 不动。
#  判（4 条；每条一行「抓什么」）：
#    R0 OpenMC 自报 SE vs 10 个种子的散布（逐层比值的中位数应在 0.8–1.25）     抓：OpenMC 那边 SE 本身不可信，重标就没意义
#    R1 ① SA1010 vs OpenMC 种子 1：4 个量 Q 落在新门槛里                       抓：开 S 后真有形状偏差
#    R2 ③ SA1011 vs OpenMC 种子 1：同上                                        抓：同上（换一个种子）
#    R3 你自己比 SA1010 vs SA1011：Q 落在新门槛里                               抓：新门槛还是不对（程序自己跟自己比都不过）
#  另外量一下 10 个种子的层间平均相关 ρ̄（只看，和 10/24 无 S 的 0.29–0.35 对照）。
#  盲区：45 对不独立（每个种子出现在 9 对里），门槛本身有涨落；以后有时间多跑种子再细标。
# ============================================================
import json
import math
import os
import itertools

from lesson11_check3 import TABLE, OUT_DIR, NL, mine_layers, shape_test

OMC = ["openmc_water_sab_s%d.json" % s for s in range(1, 11)]
MINE = {"① SA1010": "water_SA_N100000_s1010.json", "③ SA1011": "water_SA_N100000_s1011.json"}
ITEMS = TABLE[:4]          # 通量、H1 / O16 / B10 heating（(n,γ)、(n,α) 和 B10 同形状，不重复）


def omc_layers(o, t30, coe):
    return [x * coe for x in o[t30]["mean"]], [x * coe for x in o[t30]["se"]]


if __name__ == "__main__":
    for name in OMC + list(MINE.values()):
        assert os.path.exists(os.path.join(OUT_DIR, name)), "找不到 outputs\\%s" % name
    omcs = []
    for name in OMC:
        with open(os.path.join(OUT_DIR, name), "r") as f:
            o = json.load(f)
        # [Claude assert] 抓：混进了没开 S(α,β) 的、或粒子数不对的
        assert o.get("sab") == "c_H_in_H2O" and o["batches"] == 100 and o["particles"] == 100000, "%s 不是开 S、100 × 10 万" % name
        omcs.append(o)
    mine = {}
    for tag, name in MINE.items():
        with open(os.path.join(OUT_DIR, name), "r") as f:
            mine[tag] = json.load(f)

    all_ok = True
    print("%-12s %8s %8s %8s | %-22s %-22s %-22s | %s" % ("量", "R0 SE比", "新下限", "新上限", "R1 ① SA1010", "R2 ③ SA1011", "R3 自比", "ρ̄ 量出"))
    for name, sk, qk, sub, t30, ttot, mtk, mtq, coe, rho in ITEMS:
        L = [omc_layers(o, t30, coe) for o in omcs]
        n = len(L)
        # R0：自报 SE vs 种子间散布
        ratios = []
        for k in range(NL):
            mk = sum(l[0][k] for l in L) / n
            sd = math.sqrt(sum((l[0][k] - mk) ** 2 for l in L) / (n - 1))
            se = sum(l[1][k] for l in L) / n
            if sd > 0:
                ratios.append(se / sd)
        ratios.sort()
        med = ratios[len(ratios) // 2]
        ok0 = 0.8 <= med <= 1.25
        # 层间平均相关 ρ̄（只看）
        devs = []
        for l in L:
            devs.append([(l[0][k] - sum(x[0][k] for x in L) / n) / max(1e-300, math.sqrt(sum((x[0][k] - sum(y[0][k] for y in L) / n) ** 2 for x in L) / (n - 1))) for k in range(NL)])
        cs = []
        for a, b in itertools.combinations(range(NL), 2):
            cs.append(sum(d[a] * d[b] for d in devs) / (n - 1))
        rho_meas = sum(cs) / len(cs)
        # 新门槛：45 对自比
        Qs = [shape_test(L[i][0], L[i][1], L[j][0], L[j][1], rho)["Q"] for i, j in itertools.combinations(range(n), 2)]
        lo, hi = min(Qs), max(Qs)
        # R1、R2：你的 vs OpenMC 种子 1
        cells = []
        oks = [ok0]
        for tag in MINE:
            m1, s1 = mine_layers(mine[tag], sk, qk, sub)
            Q = shape_test(m1, s1, L[0][0], L[0][1], rho)["Q"]
            ok = lo <= Q <= hi
            oks.append(ok)
            cells.append("Q %6.2f %s" % (Q, "过" if ok else "✗"))
        # R3：你自己比
        a1, b1 = mine_layers(mine["① SA1010"], sk, qk, sub)
        a2, b2 = mine_layers(mine["③ SA1011"], sk, qk, sub)
        Q3 = shape_test(a1, b1, a2, b2, rho)["Q"]
        ok3 = lo <= Q3 <= hi
        oks.append(ok3)
        cells.append("Q %6.2f %s" % (Q3, "过" if ok3 else "✗"))
        all_ok = all_ok and all(oks)
        print("%-12s %5.2f %s %8.2f %8.2f | %-22s %-22s %-22s | %.3f（用 %.1f）"
              % (name, med, "过" if ok0 else "✗", lo, hi, cells[0], cells[1], cells[2], rho_meas, rho))
        print("%12s 45 对 Q 从小到大：%s" % ("", " ".join("%.0f" % q for q in sorted(Qs))))
    print("\n判据 3 丙（开 S，门槛按 OpenMC 自比 45 对重标）：%s" % ("全过 → M3 水模体对上 OpenMC" if all_ok else "有 ✗ —— 发我"))
    print("⚠️ 门槛是看了 ①③ 的数据之后才决定重标的（理由见文件头），写进 M3 小结和技术报告时照实写。")
