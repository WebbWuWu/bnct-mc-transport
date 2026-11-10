# -*- coding: utf-8 -*-
# ============================================================
#  check_tissue_angle.py  --  M3 收尾：软组织开真实角分布（H1、O16）以后，浅处 4–8% 还在不在；C12 要不要接
#  判据归 Claude，2026-11-08 写（跑之前写下门槛）。运行：Windows Ctrl+F5（只读 json，几秒钟）。
#  先跑完 explore\m3_tissue_angle_runs.py。函数、常数、OpenMC 文件表都从 check_tissue.py 拿（11/6 那份，不改它）。
#
#  比的三列（都是「每个源中子」，量和单位同 check_tissue.py 文件头）：
#    A 关 = tissue_kerma_N100000_s1010 / s1200（11/6，角分布全各向同性）
#    A 开 = tissue_kermaA_N100000_s1010 / s1200（今天，H1、O16 真实；C12、N14 仍各向同性）
#    OpenMC = M1 的 10 个种子（真实角分布、不开 S(α,β)；C12 也是真实角分布）
#
#  门槛（今天出数之前定，不因看到数据再改）
#    两段：浅处 0–3 cm = 第 1–6 层（下标 0–5）· 深处 10–15 cm = 第 21–30 层（下标 20–29），各段层加起来比
#    r = 你 / OpenMC；z = 差 / √(SE你² + SE_OMC²)；SE你 用保守的「层完全相关」（每个种子 Σse，两个种子平均再 /√2 合）
#    判：|z| < 3 → 没看到偏差 · |z| ≥ 3 且 |r − 1| > 0.5% → 偏差还在 · 其余 → 显著但 ≤ 0.5%，记下不追
#    ⚠️ 检验有没有分辨力：A 关那一列浅处如果本来就 |z| < 3，说明这个保守 SE 下看不出 4–8%，
#       「A 开 没看到偏差」不能当证据 —— 这时只看 r 往 1 靠了没有，结论写「方向对、分辨不出」
#    结论规则：
#      H1、O16 浅处 A 开都「没看到偏差」 且 C12「偏差还在」 → 建议 C12、N14 也接角分布（elastic_angle.NUCS 加两个名字 + 判据 A 重跑）
#      H1、O16、C12 浅处 A 开都「没看到偏差」              → C12、N14 不接，附录 B 记一条近似（C12、N14 质心系各向同性）
#      H1 或 O16 浅处 A 开仍「偏差还在」                  → 浅处偏高不是（只是）角分布，另查；先不接 C12
#  只看（不是判据）：两个种子各自的 r（2 个种子的散布只有 1 个自由度，不拿来算 z —— 11/7 定的「种子散布当误差」
#  要 ≥ 10 个种子才可信，这里只有 2 个，所以用保守 SE）；A 开 / A 关 逐段比；通量、D_B（角分布对它们影响应很小）
# ============================================================
import json
import math
import os

from check_tissue import TABLE, OMC_FILES, OUT_DIR, NL, DZ, N_HIST, SEEDS, Z_GATE, R_GATE, mine_one, omc_ref

BLOCKS = [("浅处 0–3 cm", list(range(0, 6))), ("深处 10–15 cm", list(range(20, 30)))]
RUNS = [("A 关", "tissue_kerma_N%d_s%d.json"), ("A 开", "tissue_kermaA_N%d_s%d.json")]


def verdict(r, z):
    if abs(z) < Z_GATE:
        return "没看到偏差"
    return "偏差还在" if abs(r - 1) > R_GATE else "显著但≤0.5%"


def load_runs(pattern):
    out = []
    for seed in SEEDS:
        p = os.path.join(OUT_DIR, pattern % (N_HIST, seed))
        assert os.path.exists(p), "找不到 %s —— 先跑 explore\\m3_tissue_angle_runs.py（A 关那份是 11/6 的 m1_tissue_runs.py）" % p
        with open(p, "r") as f:
            out.append(json.load(f))
    return out


if __name__ == "__main__":
    runs = {lab: load_runs(pat) for lab, pat in RUNS}
    omcs = []
    for name in OMC_FILES:
        with open(os.path.join(OUT_DIR, name), "r") as f:
            o = json.load(f)
        # [Claude assert] 抓：混进了开 S(α,β) 的、或粒子数不一样的 OpenMC 结果
        assert o["batches"] == 100 and o["particles"] == 10000, "%s 不是 100 批 × 1 万" % name
        omcs.append(o)
    print("你的：种子 %s（各 %d 条）× A 关 / A 开   vs   OpenMC M1 %d 个种子\n"
          % ("、".join(str(s) for s in SEEDS), N_HIST, len(omcs)))

    res = {}   # (量, 段, 列) -> (r, z, 判)
    for name, sk, qk, sub, tal, coe in TABLE:
        _, _, block = omc_ref(omcs, tal, coe)
        print("【%s】" % name)
        for bname, layers in BLOCKS:
            b, sb = block(layers)
            line = []
            sums = {}
            for lab, _ in RUNS:
                each = [mine_one(d, sk, qk, sub) for d in runs[lab]]
                a_s = [sum(e[0][k] for k in layers) for e in each]          # 每个种子这段的和
                se_s = [sum(e[1][k] for k in layers) for e in each]         # 每个种子保守 SE（层完全相关）
                a = sum(a_s) / len(a_s)
                se = math.sqrt(sum(x ** 2 for x in se_s)) / len(se_s)
                r = a / b
                z = (a - b) / math.sqrt(se ** 2 + sb ** 2)
                v = verdict(r, z)
                res[(name, bname, lab)] = (r, z, v)
                sums[lab] = a
                line.append("%s r = %.4f（种子各 %s）z = %+5.2f → %s"
                            % (lab, r, " / ".join("%.3f" % (x / b) for x in a_s), z, v))
            print("  %-12s %s" % (bname, line[0]))
            print("  %-12s %s" % ("", line[1]))
            print("  %-12s 只看：A 开 / A 关 = %.4f" % ("", sums["A 开"] / sums["A 关"]))
        print()

    print("#" * 60)
    print("结论（门槛见文件头）")
    sh = "浅处 0–3 cm"
    for name in ("H1 kerma", "O16 kerma", "C12 kerma"):
        off, on = res[(name, sh, "A 关")], res[(name, sh, "A 开")]
        print("  %-10s 浅处：A 关 r = %.4f z = %+5.2f（%s）→ A 开 r = %.4f z = %+5.2f（%s）"
              % (name, off[0], off[1], off[2], on[0], on[1], on[2]))
    no_power = [n for n in ("H1 kerma", "O16 kerma", "C12 kerma") if abs(res[(n, sh, "A 关")][1]) < Z_GATE]
    if no_power:
        print("  ⚠️ 分辨力不够：%s 在 A 关时浅处就 |z| < 3 —— 这几项「A 开没看到偏差」不算证据，只看 r 往 1 靠没靠" % "、".join(no_power))
    ok = lambda n: res[(n, sh, "A 开")][2] == "没看到偏差"
    if not (ok("H1 kerma") and ok("O16 kerma")):
        print("  → H1 或 O16 浅处开了角分布仍偏：不是（只是）角分布，另查；先不接 C12")
    elif ok("C12 kerma"):
        print("  → H1、O16、C12 浅处都没看到偏差：C12、N14 不接，附录 B 记「C12、N14 质心系各向同性」")
    else:
        print("  → H1、O16 好了、C12 还偏：建议 C12、N14 也接角分布（NUCS 加两个名字 + 判据 A 重跑）")
