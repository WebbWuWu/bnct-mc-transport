# -*- coding: utf-8 -*-
# ============================================================
#  check3_m3_seeds.py  --  M3 判据 3：用你自己的 10 个种子标定（11/6 用户定甲）
#  判据归 Claude，2026-11-06 18:00 写 —— 在 10 个种子跑出来之前写下门槛。运行：Windows Ctrl+F5（只读 json）。
#
#  为什么（照实写）：开 S 后，你的两个 10 万条种子（1010、1011）自己比，通量 Q = 70、B10 Q = 76，
#    比 OpenMC 自比 45 对的最大值（42、50）还大。猜测：你用隐式俘获 + 轮盘赌（OpenMC 默认模拟俘获），
#    一条带权重的热中子撞几百次、在好几厘米里来回走，层间相关比 OpenMC 的更强、拖得更远 → 拿 OpenMC 标的门槛判你还是偏紧。
#    要么是这个（只是相关结构不同，平均值没问题），要么是你开 S 后 SE 本身算小了（真问题）。下面 S1 把两者分开。
#
#  判（4 条；每条一行「抓什么」）
#    S1 你自报的 SE vs 10 个种子之间的散布：逐层比值的中位数在 [0.8, 1.25]（10/31 M2 30 个种子是 0.91–0.97）
#       抓：开 S 后你的 SE 算小了 / 算大了（比如历史之间不独立、计数放错层）
#    S2 门槛：你自己 10 个种子两两比 45 对的 Q → [最小, 最大] = 新门槛（只看，不判；和 OpenMC 自比 45 对并排）
#    S3 全部你的数据合起来（10 × 2 万 + 1010、1011 两个 10 万 = 40 万条）vs OpenMC 开 S 10 个种子合起来（1 亿）：
#       形状 Q 落在 S2 的门槛里                                             抓：开 S 后形状真的不一样
#    S4 同上，4 个总量 |z| ≤ 3（OpenMC 那边 SE 用 10 个种子的散布）         抓：整体偏
#  另外：逐层最大 |z|（只看）· 你和 OpenMC 各自的层间平均相关 ρ̄（只看）
#  Q 的公式、ρ̄ 都和 lesson11_check3 一样（同一个函数 shape_test），只换门槛来源。
# ============================================================
import json
import math
import os
import itertools

from lesson11_check3 import TABLE, OUT_DIR, NL, shape_test, pick

SEEDS = range(3001, 3011)
N_SMALL = 20000
BIG = {"water_SA_N100000_s1010.json": 100000, "water_SA_N100000_s1011.json": 100000}
OMC = ["openmc_water_sab_s%d.json" % s for s in range(1, 11)]
ITEMS = TABLE[:4]
TOT = {"通量": ("track_tot_sum", "track_tot_sq", None, "water total flux"),
       "H1 heating": ("kerma_tot_sum", "kerma_tot_sq", "H1", "water total H1 heating"),
       "O16 heating": ("kerma_tot_sum", "kerma_tot_sq", "O16", "water total O16 heating"),
       "B10 heating": ("kerma_tot_sum", "kerma_tot_sq", "B10", "water total B10 heating")}


def ms(s, q, n):
    """Σx、Σx²、n → 均值、均值标准误（同 mcstat.mean_se）。"""
    m = s / n
    var = n / (n - 1) * (q / n - m * m)
    return m, math.sqrt(max(var, 0.0) / n)


def layers(d, sk, qk, sub, n, coe=1.0):
    su, sq = pick(d, sk, sub), pick(d, qk, sub)
    out = [ms(su[k], sq[k], n) for k in range(NL)]
    return [x[0] * coe for x in out], [x[1] * coe for x in out]


def rho_bar(series):
    n = len(series)
    mean = [sum(s[k] for s in series) / n for k in range(NL)]
    sd = [math.sqrt(sum((s[k] - mean[k]) ** 2 for s in series) / (n - 1)) for k in range(NL)]
    cs = []
    for a, b in itertools.combinations(range(NL), 2):
        if sd[a] > 0 and sd[b] > 0:
            cs.append(sum((s[a] - mean[a]) * (s[b] - mean[b]) for s in series) / (n - 1) / sd[a] / sd[b])
    return sum(cs) / len(cs)


if __name__ == "__main__":
    small = []
    for s in SEEDS:
        p = os.path.join(OUT_DIR, "water_SA_N%d_s%d.json" % (N_SMALL, s))
        assert os.path.exists(p), "找不到 %s —— 先跑 explore\\m3_sa_seeds.py" % p
        with open(p, "r") as f:
            small.append(json.load(f))
    big = {}
    for name in BIG:
        with open(os.path.join(OUT_DIR, name), "r") as f:
            big[name] = json.load(f)
    omcs = []
    for name in OMC:
        with open(os.path.join(OUT_DIR, name), "r") as f:
            o = json.load(f)
        assert o.get("sab") == "c_H_in_H2O", "%s 没开 S" % name
        omcs.append(o)

    all_ok = True
    print("%-12s | %-14s | %-26s | %-18s | %-34s | %-10s | %s" % ("量", "S1 SE比(中位)", "S2 你自比 45 对 Q [小, 大]", "OpenMC 自比 [小, 大]", "S3 合起来 vs OpenMC 形状 Q", "S4 总量 z", "ρ̄ 你 / OpenMC · 最大|z|"))
    for name, sk, qk, sub, t30, ttot, mtk, mtq, coe, rho in ITEMS:
        L = [layers(d, sk, qk, sub, N_SMALL, coe) for d in small]
        # S1
        ratios = []
        for k in range(NL):
            mk = sum(l[0][k] for l in L) / len(L)
            sd = math.sqrt(sum((l[0][k] - mk) ** 2 for l in L) / (len(L) - 1))
            se = sum(l[1][k] for l in L) / len(L)
            if sd > 0 and mk > 0 and se / mk <= 0.2:
                ratios.append(se / sd)
        ratios.sort()
        med = ratios[len(ratios) // 2]
        ok1 = 0.8 <= med <= 1.25
        # S2
        Qs = [shape_test(L[i][0], L[i][1], L[j][0], L[j][1], rho)["Q"] for i, j in itertools.combinations(range(len(L)), 2)]
        lo, hi = min(Qs), max(Qs)
        O = [([x * coe for x in o[t30]["mean"]], [x * coe for x in o[t30]["se"]]) for o in omcs]
        Qo = [shape_test(O[i][0], O[i][1], O[j][0], O[j][1], rho)["Q"] for i, j in itertools.combinations(range(len(O)), 2)]
        # S3：你的全部合起来（Σx、Σx² 直接加）
        su = [0.0] * NL
        sq = [0.0] * NL
        n_all = 0
        for d in small:
            for k in range(NL):
                su[k] += pick(d, sk, sub)[k]
                sq[k] += pick(d, qk, sub)[k]
            n_all += N_SMALL
        for nm, d in big.items():
            for k in range(NL):
                su[k] += pick(d, sk, sub)[k]
                sq[k] += pick(d, qk, sub)[k]
            n_all += BIG[nm]
        m1 = [ms(su[k], sq[k], n_all)[0] * coe for k in range(NL)]
        s1 = [ms(su[k], sq[k], n_all)[1] * coe for k in range(NL)]
        nO = len(O)
        m2 = [sum(o[0][k] for o in O) / nO for k in range(NL)]
        s2 = [math.sqrt(sum((o[0][k] - m2[k]) ** 2 for o in O) / (nO - 1) / nO) for k in range(NL)]
        sh = shape_test(m1, s1, m2, s2, rho)
        ok3 = lo <= sh["Q"] <= hi
        zmax = max(abs(z) for z in sh["z"] if z is not None)
        # S4：总量
        tk, tq, tsub, tot_t = TOT[name]
        ts = sum(pick(d, tk, tsub) for d in small) + sum(pick(d, tk, tsub) for d in big.values())
        tsq = sum(pick(d, tq, tsub) for d in small) + sum(pick(d, tq, tsub) for d in big.values())
        tm1, tse1 = ms(ts, tsq, n_all)
        ot = [o[tot_t]["mean"][0] for o in omcs]
        tm2 = sum(ot) / nO
        tse2 = math.sqrt(sum((x - tm2) ** 2 for x in ot) / (nO - 1) / nO)
        zt = (tm1 - tm2) / math.sqrt(tse1 ** 2 + tse2 ** 2)
        ok4 = abs(zt) <= 3
        all_ok = all_ok and ok1 and ok3 and ok4
        print("%-12s | %5.2f %-8s | [%6.2f, %6.2f]           | [%6.2f, %6.2f]   | Q %6.2f %s（你/OMC 总量 %.4f）      | %+5.2f %s | %.3f / %.3f · %.2f"
              % (name, med, "过" if ok1 else "✗", lo, hi, min(Qo), max(Qo), sh["Q"], "过" if ok3 else "✗", tm1 / tm2,
                 zt, "过" if ok4 else "✗", rho_bar([l[0] for l in L]), rho_bar([o[0] for o in O]), zmax))
        print("%12s   你自比 45 对 Q 中位数 %.1f · OpenMC 自比中位数 %.1f" % ("", sorted(Qs)[22], sorted(Qo)[22]))
    print("\n判据 3（M3，按你自己的种子标定）：%s" % ("全过 → M3 水模体对上 OpenMC" if all_ok else "有 ✗ —— 发我"))
    print("⚠️ 照实写：1010 单种子判据没过之后才改成这个判法；门槛规则在 10 个种子出数之前写下。")
