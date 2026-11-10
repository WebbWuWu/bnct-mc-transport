# -*- coding: utf-8 -*-
# ============================================================
#  check_M3_seeds.py  --  判据 3 的「按你自己的种子标定」通用版（判法、门槛和 check_M3_3_seeds.py 一字不差，只换输入）
#  判据归 Claude，2026-11-06 18:50 写（用户定甲：H1、O16 形状 ✗ 查是不是角分布）。运行：Windows Ctrl+F5（默认 A）
#  或命令行 python check_M3_seeds.py A / M2 / SA。只读 json。
#
#  三种输入（你的 × OpenMC）：
#    SA = S+A 开 10 × 2 万（3001–3010）+ 1010、1011 两个 10 万  vs  OpenMC 开 S 10 个种子    ← 和 check_M3_3_seeds.py 同一个结果（自检用）
#    A  = 只开角分布 10 × 2 万（3001–3010，explore\m3_A_seeds.py 跑）  vs  OpenMC 不开 S（openmc_water_s1–10，真实角分布）
#    M2 = 都不开 30 × 2 万（2001–2030，10/31 新发生器那 30 个）        vs  同上                ← 对照组：不开 A 时形状对不对
#  怎么读：A ✗ 而 M2 过 → 偏差跟着角分布走（A 的实现，或「CM 各向同性」本来就更接近？见下）；A 过 → 问题在 S 和 A 叠起来。
#  ⚠️ OpenMC 那边永远用真实角分布，所以 M2（你各向同性）本身带着近似 #4 —— M2 过只说明在这个精度下近似 #4 看不出来。
# ============================================================
import json
import sys
import math
import os
import itertools

from check_M2_3 import TABLE, OUT_DIR, NL, shape_test, pick

CONFIG = {
    "SA": ("water_SA_N20000_s%d.json", range(3001, 3011),
           {"water_SA_N100000_s1010.json": 100000, "water_SA_N100000_s1011.json": 100000},
           "openmc_water_sab_s%d.json", "c_H_in_H2O", "先跑 explore\\m3_sa_seeds.py"),
    "A":  ("water_A_N20000_s%d.json", range(3001, 3011), {}, "openmc_water_s%d.json", None, "先跑 explore\\m3_A_seeds.py"),
    "M2": ("water_N20000_s%d.json", range(2001, 2031), {}, "openmc_water_s%d.json", None, "10/31 那 30 个种子"),
}
MODE = sys.argv[1] if len(sys.argv) > 1 else "A"
assert MODE in CONFIG, "参数只能是 SA / A / M2，收到 %s" % MODE
PAT, SEEDS, BIG, OMC_PAT, SAB, HOW = CONFIG[MODE]
N_SMALL = 20000
OMC = [OMC_PAT % s for s in range(1, 11)]
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
        p = os.path.join(OUT_DIR, PAT % s)
        assert os.path.exists(p), "找不到 %s —— %s" % (p, HOW)
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
        assert o.get("sab") == SAB, "%s 的 sab = %s，这一组应该是 %s" % (name, o.get("sab"), SAB)
        omcs.append(o)

    all_ok = True
    print("%-12s | %-14s | %-26s | %-18s | %-34s | %-10s | %s" % ("量", "S1 SE比(中位)", "S2 你自比 %d 对 Q [小, 大]" % (len(SEEDS) * (len(SEEDS) - 1) // 2), "OpenMC 自比 [小, 大]", "S3 合起来 vs OpenMC 形状 Q", "S4 总量 z", "ρ̄ 你 / OpenMC · 最大|z|"))
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
        try:
            Qs = [shape_test(L[i][0], L[i][1], L[j][0], L[j][1], rho)["Q"] for i, j in itertools.combinations(range(len(L)), 2)]
        except KeyError:
            print("%-12s | 不判：单个 2 万条种子 R ≤ 0.2 的层不到 16 层，Q 的自由度对不齐（M2 的 O16 深处）" % name)
            continue
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
        print("%12s   你自比 Q 中位数 %.1f · OpenMC 自比中位数 %.1f" % ("", sorted(Qs)[len(Qs) // 2], sorted(Qo)[len(Qo) // 2]))
    print("\n[%s] 判据 3（按你自己的种子标定）：%s" % (MODE, "全过" if all_ok else "有 ✗ —— 发我"))
