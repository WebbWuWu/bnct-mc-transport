# -*- coding: utf-8 -*-
# ============================================================
#  lesson11_check_rho.py  --  M2 判据 3（方案丙）的前置：用 OpenMC 水模体 10 个种子量层间相关 ρ̄，并查 OpenMC 的 SE
#  判据归 Claude，2026-10-24 写。被测：lesson11_openmc.py 跑出的 outputs\openmc_water_s1–s10.json（100 批 × 10 万）。
#  运行：Windows Ctrl+F5 或 WSL 都行（只读 json，不 import openmc）。几秒钟。
#  🔴 规则 14：这个脚本【只看 OpenMC 自己】，不碰你的程序的结果 —— ρ̄ 和门槛要在两边比之前定下来。
# ============================================================
#
#  为什么要量 ρ̄（判据 3 重标方案 第 3、4 节）
#    方案丙的「形状」检验：每层 z 减掉 30 层的平均 z̄，Σ(z − z̄)² / (1 − ρ̄) 按 χ²₂₉ 判，门槛 [13.12, 52.34]（99%）。
#    ρ̄ 是层与层 z 的相关系数平均。ρ̄ 估大比估小伤得多（估大 0.1 → 都对时误报 4.1%；估小 0.1 → 1.8%），
#    所以 10 个种子量出来以后【往小取到 0.1 的整数倍】。
#
#  记号
#    M      独立运行数（10 个种子，同一设置）
#    s_k    第 k 层：10 次运行均值的样本标准差（路 2：实际散多开）
#    se_k   第 k 层：10 次运行自报 SE 的均方根（路 1：OpenMC 说自己散多开）
#    R      sqrt(30 层 (s_k / se_k)² 的平均)，SE 对的话 ≈ 1
#    ρ̄      30 层两两相关系数（435 对）的平均
#
#  判据（4 条；每条一行「抓什么」）
#    W1 10 个 json 的种子各不相同、粒子数都是 10 万 × 100 批，结果两两不完全相同
#       抓：种子没传进去（10 次一模一样 → 散布 = 0，下面全是假的）
#    W2 每个「总量」tally = 同名 30 层 tally 加起来（相对差 ≤ 1e-9；mesh 正好盖住整块板，两者应该是同一组径迹）
#       抓：不分层的 tally 漏了核素 / 名字配错、mesh 没盖满板
#    W3 每个 30 层 tally 的 R 在 [0.44, 1.62]（9 个自由度、30 层完全相关的最坏情况，99%；来历同 lesson09_check_se.py）
#       另外每个总量 tally：10 次的样本标准差 / 自报 SE 在 [sqrt(χ²₉(0.005)/9), sqrt(χ²₉(0.995)/9)] = [0.44, 1.62]
#       抓：OpenMC 的 SE 估小 / 估大一倍以上（那样 z 全是虚的）
#    W4 试响：把路 1 的 SE 故意 × 0.5，6 个 30 层 tally 都必须出界
#       抓：门槛太松
#  只量不判：每个 tally 的 ρ̄ → 定下方案丙要用的 ρ̄
# ============================================================
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "outputs")
NL = 30
SEEDS = list(range(1, 11))
LAYERED = ["water flux", "water H1 heating", "water O16 heating", "water B10 heating", "water H1 n,gamma", "water B10 n,a"]
TOTALS = {"water total flux": "water flux", "water total H1 heating": "water H1 heating",
          "water total O16 heating": "water O16 heating", "water total B10 heating": "water B10 heating"}
R_LO, R_HI = 0.44, 1.62
SHAPE_LO, SHAPE_HI = 13.12, 52.34        # χ²₂₉ 的 0.5% / 99.5% 分位（方案文档第 3 节）


def sample_sd(x):
    m = sum(x) / len(x)
    return math.sqrt(sum((v - m) ** 2 for v in x) / (len(x) - 1))


def ratio_R(runs, name, f=1.0):
    rs = []
    for k in range(NL):
        s_k = sample_sd([o[name]["mean"][k] for o in runs])
        se_k = math.sqrt(sum((f * o[name]["se"][k]) ** 2 for o in runs) / len(runs))
        rs.append(s_k / se_k)
    return math.sqrt(sum(r * r for r in rs) / NL)


def mean_corr(runs, name):
    M = len(runs)
    cols = [[o[name]["mean"][k] for o in runs] for k in range(NL)]
    mus = [sum(c) / M for c in cols]
    sds = [sample_sd(c) for c in cols]
    acc, cnt = 0.0, 0
    for i in range(NL):
        for j in range(i + 1, NL):
            cov = sum((cols[i][t] - mus[i]) * (cols[j][t] - mus[j]) for t in range(M)) / (M - 1)
            acc += cov / (sds[i] * sds[j])
            cnt += 1
    return acc / cnt


if __name__ == "__main__":
    runs = []
    for s in SEEDS:
        with open(os.path.join(OUT_DIR, "openmc_water_s%d.json" % s), "r") as f:
            runs.append(json.load(f))

    # ---- W1 ----
    seeds = [o["seed"] for o in runs]
    assert sorted(seeds) == SEEDS and all(o["particles"] == 100000 and o["batches"] == 100 for o in runs), \
        "W1：种子 %r / 粒子数、批数不对" % (seeds,)
    firsts = [o["water H1 heating"]["mean"][0] for o in runs]
    assert len(set(firsts)) == len(firsts), "W1：有两次运行第 1 层 H1 heating 完全相同 —— 种子没传进去？"
    print("W1 过：10 个种子 %s，各 100 批 × 10 万，结果两两不同" % seeds)

    # ---- W2 ----
    for tot, lay in TOTALS.items():
        for o in runs:
            a, b = o[tot]["mean"][0], sum(o[lay]["mean"])
            assert abs(a - b) <= 1e-9 * abs(b), "W2：种子 %d 的 %s = %.12g，30 层之和 %.12g" % (o["seed"], tot, a, b)
    print("W2 过：4 个总量 tally = 同名 30 层 tally 之和（10 个种子都是）")

    # ---- W3 + 量 ρ̄ ----
    print("\n%-22s %8s %8s   %s" % ("tally", "R", "ρ̄", "30 层 R=SE/均值 最大"))
    rhos = {}
    for name in LAYERED:
        R = ratio_R(runs, name)
        rho = mean_corr(runs, name)
        rhos[name] = rho
        rel = max(o[name]["se"][k] / o[name]["mean"][k] for o in runs for k in range(NL))
        print("%-22s %8.3f %8.3f   %.4f" % (name, R, rho, rel))
        assert R_LO <= R <= R_HI, "W3：%s 的 R = %.3f 不在 [%.2f, %.2f] —— OpenMC 的 SE 不可信" % (name, R, R_LO, R_HI)
    for tot in TOTALS:
        vals = [o[tot]["mean"][0] for o in runs]
        se = math.sqrt(sum(o[tot]["se"][0] ** 2 for o in runs) / len(runs))
        q = sample_sd(vals) / se
        print("%-26s 10 次散布 / 自报 SE = %.3f   （均值 %.6g，自报 R %.5f）" % (tot, q, sum(vals) / len(vals), se / (sum(vals) / len(vals))))
        assert R_LO <= q <= R_HI, "W3：%s 的 散布/SE = %.3f 不在 [%.2f, %.2f]" % (tot, q, R_LO, R_HI)
    print("W3 过：6 个 30 层 tally、4 个总量 tally 的 SE 都可信")

    # ---- W4 ----
    for name in LAYERED:
        assert not (R_LO <= ratio_R(runs, name, 0.5) <= R_HI), "W4：%s 的 SE × 0.5，W3 居然没响" % name
    print("W4 试响过：SE × 0.5 → 6 个 tally 的 R 都出界")

    # ---- 定 ρ̄ ----
    print("\n方案丙要用的 ρ̄（往小取到 0.1 的整数倍；ρ̄ < 0 按 0 算）：")
    for name in LAYERED:
        use = max(0.0, math.floor(rhos[name] * 10) / 10)
        print("  %-22s 量出 %.3f → 用 %.1f   形状门槛：Σ(z − z̄)² / %.1f ∈ [%.2f, %.2f]"
              % (name, rhos[name], use, 1 - use, SHAPE_LO, SHAPE_HI))
