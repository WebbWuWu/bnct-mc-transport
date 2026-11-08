# -*- coding: utf-8 -*-
# ============================================================
#  check_tissue_angle10.py  --  软组织开角分布（H1、O16）× 10 个种子 vs OpenMC M1：浅处那约 4% 到底在不在
#  判据归 Claude，2026-11-08 写（11/8 用户定甲：补 8 个种子；门槛在出数之前写下，不因看到数据再改）。
#  先跑完 explore\m3_tissue_angle_runs.py（10 个种子：1010、1200、2001–2008）。运行：Windows Ctrl+F5（几秒钟）。
#
#  为什么换误差的算法：check_tissue_angle.py 只有 2 个种子，种子散布只有 1 个自由度不能用，只好用「层完全相关」的
#  保守 SE —— 结果 A 关时浅处就 |z| < 3，分辨不出 4–6%。现在 10 个种子 → 照 11/7 定的规矩：
#  多层合成一个数，SE 一律用种子之间的散布（层间相关自动在里面，不用猜）。
#    你：10 个种子各算一个「段和」→ 平均、SE = 标准差 / √10 · OpenMC：同样，M1 的 10 个种子（check_tissue.omc_ref 的 block）
#
#  门槛
#    两段：浅处 0–3 cm（下标 0–5）· 深处 10–15 cm（下标 20–29），量 = 通量、D_B、H1 / C12 / O16 kerma
#    r = 你 / OpenMC，z = 差 / √(SE你² + SE_OMC²)，分辨力 MDD = 3 × √(SE你² + SE_OMC²) / OpenMC（能看出的最小相对偏差）
#    判：|z| ≥ 3 且 |r − 1| > 0.5% → 偏差还在
#        |z| ≥ 3 且 |r − 1| ≤ 0.5% → 显著但 ≤ 0.5%，记下不追
#        |z| < 3 且 MDD ≤ 2%      → 没看到偏差（这次分辨得出，算数）
#        |z| < 3 且 MDD > 2%      → 分辨不出（再加种子也只能写「≤ MDD」）
#    浅处结论：
#      H1、C12、O16 浅处都「没看到偏差」              → 浅处 4% 是 2 个种子的涨落，结案；C12、N14 不接角分布 = 近似 #39
#      三个都「偏差还在」且 r 差不多（彼此差 < 2%）     → 三者共有的偏高 = 中子能谱 / 源的问题，不是角分布；记已知问题另查
#      只有 C12「偏差还在」                          → 回头考虑给 C12 接角分布
#  只看：每个种子各自的 r（看有没有离群的种子）
# ============================================================
import json
import math
import os

from check_tissue import TABLE, OMC_FILES, OUT_DIR, N_HIST, Z_GATE, R_GATE, mine_one, omc_ref

SEEDS10 = (1010, 1200) + tuple(range(2001, 2009))
BLOCKS = [("浅处 0–3 cm", list(range(0, 6))), ("深处 10–15 cm", list(range(20, 30)))]
MDD_GATE = 0.02


def verdict(r, z, mdd):
    if abs(z) >= Z_GATE:
        return "偏差还在" if abs(r - 1) > R_GATE else "显著但≤0.5%"
    return "没看到偏差" if mdd <= MDD_GATE else "分辨不出"


if __name__ == "__main__":
    mines = []
    for seed in SEEDS10:
        p = os.path.join(OUT_DIR, "tissue_kermaA_N%d_s%d.json" % (N_HIST, seed))
        assert os.path.exists(p), "找不到 %s —— 先跑完 explore\\m3_tissue_angle_runs.py" % p
        with open(p, "r") as f:
            mines.append(json.load(f))
    omcs = []
    for name in OMC_FILES:
        with open(os.path.join(OUT_DIR, name), "r") as f:
            o = json.load(f)
        # [Claude assert] 抓：混进了开 S(α,β) 的、或粒子数不一样的 OpenMC 结果
        assert o["batches"] == 100 and o["particles"] == 10000, "%s 不是 100 批 × 1 万" % name
        omcs.append(o)
    n = len(mines)
    print("你的：开角分布（H1、O16）%d 个种子 × %d 条   vs   OpenMC M1 %d 个种子；SE 两边都用种子散布\n"
          % (n, N_HIST, len(omcs)))

    res = {}
    for name, sk, qk, sub, tal, coe in TABLE:
        _, _, block = omc_ref(omcs, tal, coe)
        print("【%s】" % name)
        for bname, layers in BLOCKS:
            b, sb = block(layers)
            sums = [sum(mine_one(d, sk, qk, sub)[0][k] for k in layers) for d in mines]
            a = sum(sums) / n
            sa = math.sqrt(sum((x - a) ** 2 for x in sums) / (n - 1) / n)
            r = a / b
            se = math.sqrt(sa ** 2 + sb ** 2)
            z = (a - b) / se
            mdd = 3 * se / b
            v = verdict(r, z, mdd)
            res[(name, bname)] = (r, z, mdd, v)
            print("  %-12s r = %.4f  z = %+5.2f  分辨力 %.1f%%  → %s" % (bname, r, z, 100 * mdd, v))
            print("  %-12s 只看，各种子 r：%s" % ("", " ".join("%.3f" % (x / b) for x in sums)))
        print()

    print("#" * 60)
    sh = "浅处 0–3 cm"
    ks = ("H1 kerma", "C12 kerma", "O16 kerma")
    for k in ks:
        r, z, mdd, v = res[(k, sh)]
        print("  %-10s 浅处 r = %.4f z = %+5.2f 分辨力 %.1f%% → %s" % (k, r, z, 100 * mdd, v))
    vs = [res[(k, sh)][3] for k in ks]
    rs = [res[(k, sh)][0] for k in ks]
    if all(v == "没看到偏差" for v in vs):
        print("  → 浅处 4% 是 2 个种子的涨落，结案；C12、N14 不接角分布 = 近似 #39")
    elif all(v == "偏差还在" for v in vs) and max(rs) - min(rs) < 0.02:
        print("  → 三个核素一起偏高约 %.1f%%：能谱 / 源的问题，不是角分布；记已知问题另查，C12 不接" % (100 * (sum(rs) / 3 - 1)))
    elif vs[1] == "偏差还在" and vs[0] != "偏差还在" and vs[2] != "偏差还在":
        print("  → 只有 C12 还偏：回头考虑给 C12 接角分布")
    else:
        print("  → 不在预先写好的三种情况里：把整段输出发 Claude，照实记，不临时改门槛")
