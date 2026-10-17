# -*- coding: utf-8 -*-
# ============================================================
#  lesson09_check_se.py  --  第 9 课判据 4 的前置：OpenMC 的 SE（误差棒）可信不可信
#  判据归 Claude，2026-10-17 写。被测：OpenMC 自由气体 KUR 真问题，100 批 × 1 万，7 个 tally。
#  运行：🔴 只在 WSL 里手动跑（要 import openmc）：
#        先输 wsl → conda activate /home/siyiwu/miniconda3/envs/openmc
#        → cd /mnt/e/BNCT_GPU/my_code → python lesson09_check_se.py
#        第一次会补跑 8 个种子（5–12），每个约 10 s；已经有 json 的种子直接读、不重跑。
# ============================================================
#
#  为什么要查（10/16 判据 5 的 5-2）
#    S(α,β) 种子 3 vs 4 互比：O16 heating χ² = 74、H1 48、C12 51（30 层，应 ≈ 30）。
#    判据 4 比的就是 heating —— 要是 OpenMC 的 SE 估小了，判据 4 的 z 就全是虚的。
#
#  两个假说
#    (A) SE 是对的，只是 30 层的涨落是「一起动」的（相邻层强相关）：
#        快中子少数几个高能的决定了大部分 heating，某一次运行源里高能快中子多抽了几个，
#        30 层就一起偏高。这时 χ² 的【平均】仍是 30，但单次 χ² 会散得很开（74 不稀奇）。
#    (B) SE 本身估小了（批均值分布拖尾等）—— 这时 χ² 的【平均】就会明显大于 30。
#    10/17 开工前用现有 4 份 json 预查：种子 3 的 H1 heating 30 层总和比其他三次低 ~1%，
#    种子 3 一出场，30 层 z 就 27/30 同号 —— 像 (A)。但只有 4 次运行，不够下结论 → 补 8 个种子。
#
#  两条独立的路（纪律 2）
#    路 1：OpenMC 每次运行里自己报的 SE（来自 100 批之间的散布）
#    路 2：10 次独立运行（种子不同）的结果彼此之间实际散多开（样本标准差）
#    比的是同一个量：「一次 100 批 × 1 万的运行，结果的标准差」
#
#  记号
#    M        独立运行次数（= 10：种子 1、2、5–12，全是自由气体）
#    s_k      第 k 层：M 次运行均值的样本标准差（路 2）
#    se_k     第 k 层：M 次运行自己报的 SE 的均方根（路 1）
#    r_k      s_k / se_k，SE 对的话 ≈ 1
#    R        sqrt(30 层 r_k² 的平均)：一个 tally 一个数
#    ρ̄        层与层之间相关系数的平均（M 次运行算的，只看量级）
#
#  判据（≤ 5 条，每条旁边写抓什么）
#    SE-1 这 10 次真是「同一个设置、10 个不同种子」：
#         materials.xml 里都没有 <sab；settings.xml 的种子 10 个各不相同；结果两两不完全相同
#         抓：补跑时种子没传进去（10 次结果一模一样 → 散布 = 0 → 路 2 是假的）、混进了 S(α,β)
#    SE-2 每个 tally 的 R 在 [0.44, 1.62] 里
#         门槛来历：最坏情况 30 层完全相关，只剩 M − 1 = 9 个自由度，
#         sqrt(χ²₉ / 9) 的 99% 区间 = [sqrt(1.735/9), sqrt(23.59/9)] = [0.44, 1.62]
#         抓：SE 估小（或估大）一倍以上 —— 假说 (B)
#    SE-3 判据自己会响（纪律 4）：把路 1 的 SE 故意 × 0.5 再算一遍 R，7 个 tally 都必须出界
#         抓：门槛太松、什么都放过
#    SE-4 只量、不判：ρ̄、30 层总和的散布、45 对两两互比的 χ²（平均 / 最小 / 最大 / 出界几对）
#         —— 假说 (A) 的证据；也决定判据 4 能不能照搬判据 3 的「χ² ∈ [13.79, 53.67]」
# ============================================================

import math
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "outputs")
RUN_DIR = os.path.join(OUT_DIR, "openmc")

NL = 30
OLD = [(1, "openmc_kur.json", "kur"), (2, "openmc_kur_s2.json", "kur_s2")]
NEW_SEEDS = [5, 6, 7, 8, 9, 10, 11, 12]
RUNS = OLD + [(s, "openmc_kur_s%d.json" % s, "kur_s%d" % s) for s in NEW_SEEDS]

NAMES = [("flux", "kur flux"), ("D_B", "kur B10 n,a"), ("D_N", "kur N14 n,p"), ("(n,g)", "kur n,gamma"),
         ("H1 heat", "kur H1 heating"), ("C12 heat", "kur C12 heating"), ("O16 heat", "kur O16 heating")]
R_LO, R_HI = 0.44, 1.62
CHI2_LO, CHI2_HI = 13.79, 53.67      # 判据 3 用的 30 自由度 99% 区间


def run_missing():
    """补跑还没有 json 的种子（只在装了 openmc 的环境里）。"""
    todo = [r for r in RUNS if not os.path.exists(os.path.join(OUT_DIR, r[1]))]
    if not todo:
        print("  10 份 json 都在，不重跑")
        return
    try:
        import openmc  # noqa: F401
    except ImportError:
        raise SystemExit("还缺 %d 份 json（种子 %s），但这个环境没有 openmc —— 请在 WSL 里跑"
                         % (len(todo), [r[0] for r in todo]))
    from lesson09_openmc import kur, tru_tally, mater, run, tallies, KUR_TALLIES
    for seed, fname, d in todo:
        print("  补跑种子 %d → %s" % (seed, fname))
        sp = run(kur(), tru_tally(), os.path.join(RUN_DIR, d), batch=100, n_bat=10000, mat=mater(), seed=seed)
        out = {}
        for a, b, c in KUR_TALLIES:
            m, s = tallies(sp, a)
            out[a] = {"mean": m, "se": s}
        out["batches"] = 100
        out["particles"] = 10000
        out["seed"] = seed
        with open(os.path.join(OUT_DIR, fname), "w") as f:
            json.dump(out, f)


def seed_in_xml(d):
    p = os.path.join(RUN_DIR, d, "settings.xml")
    assert os.path.exists(p), "找不到 %s" % p
    with open(p, "r") as f:
        txt = f.read()
    m = re.search(r"<seed>(\d+)</seed>", txt)
    return int(m.group(1)) if m else 1          # 不写 <seed> 时 OpenMC 默认种子 1


def sample_sd(xs):
    n = len(xs)
    mu = sum(xs) / n
    return math.sqrt(sum((x - mu) ** 2 for x in xs) / (n - 1))


def ratio_R(runs, name, se_factor=1.0):
    """返回 (R, 30 个 r_k)。se_factor 只给 SE-3 故意写错用。"""
    rs = []
    for k in range(NL):
        s_k = sample_sd([o[name]["mean"][k] for o in runs])
        se_k = math.sqrt(sum((se_factor * o[name]["se"][k]) ** 2 for o in runs) / len(runs))
        rs.append(s_k / se_k)
    return math.sqrt(sum(r * r for r in rs) / NL), rs


def mean_corr(runs, name):
    """层与层相关系数的平均（M 次运行算，只看量级）。"""
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


def pair_chi2(a, b, name):
    return sum((a[name]["mean"][k] - b[name]["mean"][k]) ** 2 /
               (a[name]["se"][k] ** 2 + b[name]["se"][k] ** 2) for k in range(NL))


if __name__ == "__main__":
    print("=" * 78)
    print("准备：补跑缺的种子")
    run_missing()
    runs = []
    for seed, fname, d in RUNS:
        with open(os.path.join(OUT_DIR, fname), "r") as f:
            runs.append(json.load(f))

    # ---------------- SE-1 ----------------
    print("=" * 78)
    print("SE-1 10 次是同一个设置、10 个不同种子")
    seeds = []
    for (seed, fname, d), o in zip(RUNS, runs):
        p = os.path.join(RUN_DIR, d, "materials.xml")
        assert os.path.exists(p), "找不到 %s" % p
        with open(p, "r") as f:
            assert "<sab" not in f.read(), "%s 里有 S(α,β)，这一组只该有自由气体" % p
        assert o.get("batches") == 100 and o.get("particles") == 10000, "%s 不是 100 批 × 1 万" % fname
        seeds.append(seed_in_xml(d))
    print("  OpenMC 实际用的种子：%s" % seeds)
    assert len(set(seeds)) == len(RUNS), "种子有重复：%s —— 补跑时 seed 没传进 run()？" % seeds
    for i in range(len(runs)):
        for j in range(i + 1, len(runs)):
            assert runs[i]["kur flux"]["mean"] != runs[j]["kur flux"]["mean"], \
                "第 %d、%d 次（种子 %d、%d）通量逐层完全相同 —— 其实是同一次运行" % (i, j, seeds[i], seeds[j])
    print("  SE-1 过")

    # ---------------- SE-2 ----------------
    print("=" * 78)
    print("SE-2 路 2（10 次之间实际散布）/ 路 1（OpenMC 自报 SE）：R 应在 [%.2f, %.2f]" % (R_LO, R_HI))
    print("  %-10s %6s   r_k 的范围        第 1、6、11、16、21、26 层的 r_k" % ("量", "R"))
    Rs = {}
    for disp, name in NAMES:
        R, rs = ratio_R(runs, name)
        Rs[name] = R
        print("  %-10s %6.2f   [%.2f, %.2f]   %s" % (disp, R, min(rs), max(rs), " ".join("%.2f" % r for r in rs[::5])))
    for disp, name in NAMES:
        assert R_LO <= Rs[name] <= R_HI, \
            "%s：10 次之间实际散布是 OpenMC 自报 SE 的 %.2f 倍（应在 [%.2f, %.2f]）—— SE 估%s了" \
            % (disp, Rs[name], R_LO, R_HI, "小" if Rs[name] > R_HI else "大")
    print("  SE-2 过：OpenMC 自报的 SE 和实际散布对得上")

    # ---------------- SE-3 ----------------
    print("=" * 78)
    print("SE-3 判据自己会响：把 SE 故意 × 0.5（相当于 OpenMC 把 SE 估小一半）")
    for disp, name in NAMES:
        R_bad, _ = ratio_R(runs, name, se_factor=0.5)
        print("  %-10s R = %.2f" % (disp, R_bad))
        assert not (R_LO <= R_bad <= R_HI), "%s：SE 小一半时 R = %.2f 仍在门槛内 —— SE-2 太松" % (disp, R_bad)
    print("  SE-3 过：SE 小一半时 7 个都响")

    # ---------------- SE-4 ----------------
    print("=" * 78)
    print("SE-4 只量：层间相关、30 层总和的散布、45 对两两互比的 χ²")
    print("  %-10s %6s  %22s  %32s" % ("量", "ρ̄", "总和散布 / 自报（层独立 ~ 全相关）", "45 对 χ²：平均 [最小, 最大] 出界几对"))
    for disp, name in NAMES:
        rho = mean_corr(runs, name)
        tots = [sum(o[name]["mean"]) for o in runs]
        mu = sum(tots) / len(tots)
        sd_tot = sample_sd(tots) / mu
        se_ind = sum(math.sqrt(sum(s * s for s in o[name]["se"])) for o in runs) / len(runs) / mu
        se_full = sum(sum(o[name]["se"]) for o in runs) / len(runs) / mu
        chis = [pair_chi2(runs[i], runs[j], name) for i in range(len(runs)) for j in range(i + 1, len(runs))]
        n_out = sum(1 for c in chis if not (CHI2_LO <= c <= CHI2_HI))
        print("  %-10s %+6.2f  %6.2f%%  (%.2f%% ~ %.2f%%)  %10.1f [%5.1f, %5.1f]  %2d / %d"
              % (disp, rho, 100 * sd_tot, 100 * se_ind, 100 * se_full, sum(chis) / len(chis),
                 min(chis), max(chis), n_out, len(chis)))
    print()
    print("  怎么读：")
    print("    · 45 对 χ² 的【平均】≈ 30 → SE 对（和 SE-2 是一回事，换个角度看）")
    print("    · ρ̄ 大、「出界几对」多、总和散布靠近「全相关」那一端 → 假说 (A)：30 层一起动")
    print("      → 判据 4 不能照搬「χ² ∈ [13.79, 53.67]」，门槛要另定（框架问题）")
    print("    ⚠️ 45 对共用 10 次运行，彼此不独立；ρ̄ 只用 10 个样本算，误差 ±0.3 左右，只看量级")
    print("=" * 78)
    print("SE-1、SE-2、SE-3 过 → OpenMC 的 heating 误差棒可以拿来做判据 4")
