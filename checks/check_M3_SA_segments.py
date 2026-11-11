# -*- coding: utf-8 -*-
# ============================================================
#  check_M3_SA_segments.py  --  M3：S(α,β) + 角分布都开的那组，按深度分 6 段和 OpenMC 比
#  判据归 Claude，2026-11-11 夜写。只读 json，几秒钟。运行：Windows Ctrl+F5。
#
#  起因：11/11 画 M3 图时看到 S+A 那 40 万条（3001–3010 + 1010、1011）在 10–15 cm 一段
#        H1 / O16 kerma ÷ OpenMC = 1.090 / 1.108（z = +2.3 / +2.7）。那是事后看数据看到的，没有事先定门槛，不能下结论。
#        用户 11/11 23:58 定：补 20 个种子，按下面的规则判。
#
#  🔴 判定规则（2026-11-11 写下，新种子 3011–3030 一个都还没跑；出数以后不许改）：
#     用 30 个 2 万条的种子（3001–3030）。**不用** 1010、1011 那两个 10 万条 —— 偏高就是在它们身上看到的，
#     再拿来判等于同一份数据用两遍。（3001–3010 也是看过的：它们 10–15 cm 平均 1.036，不显著；留着用，照实说明。）
#     只判 10–15 cm 这一段的 H1、O16 kerma（本程序 ÷ OpenMC 开 S），误差 = 种子之间的散布：
#        两个都 |z| ≤ 2          → 涨落，结案
#        任意一个 |z| ≥ 3        → 真问题：写成已知偏差，去查「S 和 A 同时开」才走到的代码路径
#        其余                    → 「未定」，照实记，不结案
#     其余 5 段和 B10 只报数；有 |z| > 2 的照样要报出来，但不进上面的判定。
#     H1 和 O16 由同一批快中子决定，两个 z 高度相关，只算一条证据（所以规则里是「两个都」/「任意一个」，不是数个数）。
#
#  为什么用种子散布当误差：一段里好几层是一起涨落的，各层标准误平方相加会把误差低估一半多（10/31 的教训）。
#  为什么不看形状 Q：M3 小结第六节第 3 条 —— 自标门槛的 Q 会把「太好」也判 ✗；分段 z 更直接。
#  种子不足 30 个时：只打印「还差几个」，**一个数都不打印**（免得跑到一半先偷看）。
# ============================================================
# ---- 路径：本文件在 my_code 的子文件夹里，把 my_code 根目录和各子文件夹加进搜索路径 ----
import os, sys
_MY_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # my_code 根目录
for _sub in ("whiteboard", "lessons", "openmc_compare", "checks", ""):
    _p = os.path.join(_MY_CODE, _sub) if _sub else _MY_CODE
    if _p not in sys.path:
        sys.path.insert(0, _p)

import json
import math

OUT_DIR = os.path.join(_MY_CODE, "outputs")
N_SEED = 20000
SEEDS = list(range(3001, 3031))
OMC = ["openmc_water_sab_s%d.json" % i for i in range(1, 11)]
NL = 30
SEGS = [(0, 4), (4, 8), (8, 12), (12, 16), (16, 20), (20, 30)]      # 层号 [起, 止)；每层 0.5 cm → 0–2 / 2–4 / 4–6 / 6–8 / 8–10 / 10–15 cm
JUDGE_SEG = (20, 30)
JUDGE_NUCS = ("H1", "O16")
NUCS = ("H1", "O16", "B10")


def mean_sem(xs):
    """一组数的均值、均值的标准误（样本标准差 / √n）。"""
    n = len(xs)
    m = sum(xs) / n
    sd = math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1))
    return m, sd / math.sqrt(n)


def segment_table(mine, omcs):
    """mine：每个种子的 json（字典）；omcs：OpenMC 每个种子的 json。返回 {核素: [(起, 止, 比值, z, 我的相对误差), …]}。"""
    out = {}
    for nuc in NUCS:
        rows = []
        for a, b in SEGS:
            xs = [sum(d["kerma_sum"][nuc][a:b]) / N_SEED for d in mine]          # 每个种子：这一段每条历史的平均 kerma（eV·cm 的内部单位，两边同乘一个常数，比值里消掉）
            os_ = [sum(o["water %s heating" % nuc]["mean"][a:b]) for o in omcs]
            m1, s1 = mean_sem(xs)
            m2, s2 = mean_sem(os_)
            rows.append((a, b, m1 / m2, (m1 - m2) / math.sqrt(s1 ** 2 + s2 ** 2), s1 / m1))
        out[nuc] = rows
    return out


def verdict(table):
    zs = {nuc: [r[3] for r in table[nuc] if (r[0], r[1]) == JUDGE_SEG][0] for nuc in JUDGE_NUCS}
    if all(abs(z) <= 2 for z in zs.values()):
        return "涨落，结案", zs
    if any(abs(z) >= 3 for z in zs.values()):
        return "真问题", zs
    return "未定", zs


def _selftest():
    """自检：造假数据（不含随机数，结果是确定的），确认三种判定都走得到。不读真文件。
    每个假种子的值 = 1 ± 5%（正负交替，30 个平均正好是 1），OpenMC = 1 ± 0.2%；深处 10–15 cm 再整体乘 (1 + shift)。
    种子散布给出的 SE ≈ 0.05 × √(30/29) / √30 ≈ 0.93% → shift = 0 时 z = 0；2.35% 时 z ≈ 2.5；5% 时 z ≈ 5。"""
    def fake(shift_deep, n=30):
        mine = []
        for i in range(n):
            f = 1.05 if i % 2 == 0 else 0.95
            mine.append({"kerma_sum": {nuc: [N_SEED * f * (1.0 + (shift_deep if k >= 20 else 0.0)) for k in range(NL)] for nuc in NUCS}})
        omcs = [{"water %s heating" % nuc: {"mean": [1.002 if j % 2 == 0 else 0.998 for k in range(NL)]} for nuc in NUCS} for j in range(10)]
        return mine, omcs
    for shift, want in ((0.0, "涨落，结案"), (0.0235, "未定"), (0.05, "真问题")):
        got, zs = verdict(segment_table(*fake(shift)))
        assert got == want, "自检：深处偏 %+.2f%% 的假数据应判「%s」，判成了「%s」（z = %s）" % (shift * 100, want, got, zs)
    print("自检过：深处偏 0 → 结案；+2.35% → 未定；+5% → 真问题")


if __name__ == "__main__":
    _selftest()
    have = [s for s in SEEDS if os.path.exists(os.path.join(OUT_DIR, "water_SA_N%d_s%d.json" % (N_SEED, s)))]
    if len(have) < len(SEEDS):
        print("现在有 %d / %d 个种子，还差 %d 个 → 先跑 explore\\m3_sa_seeds.py（已有的自动跳过）。种子不齐不出数。" % (len(have), len(SEEDS), len(SEEDS) - len(have)))
        sys.exit(0)
    mine = []
    for s in SEEDS:
        with open(os.path.join(OUT_DIR, "water_SA_N%d_s%d.json" % (N_SEED, s)), "r") as f:
            mine.append(json.load(f))
    omcs = []
    for name in OMC:
        with open(os.path.join(OUT_DIR, name), "r") as f:
            o = json.load(f)
        assert o.get("sab") == "c_H_in_H2O", "%s 不是开 S(α,β) 的那组 OpenMC" % name
        omcs.append(o)
    table = segment_table(mine, omcs)
    print("\nS(α,β) + 角分布，%d 个种子 × %d 条 = %d 万条 vs OpenMC 开 S 10 个种子；误差 = 种子散布" % (len(SEEDS), N_SEED, len(SEEDS) * N_SEED // 10000))
    print("%-10s" % "深度 cm" + "".join("  %-22s" % (nuc + " ÷ OpenMC（z）") for nuc in NUCS))
    for i, (a, b) in enumerate(SEGS):
        line = "%4.1f–%-5.1f" % (a * 0.5, b * 0.5)
        for nuc in NUCS:
            _, _, r, z, rel = table[nuc][i]
            line += "  %.3f（%+.1f）%s        " % (r, z, " ⚠" if abs(z) > 2 else "  ")
        print(line)
    v, zs = verdict(table)
    print("\n判定（只看 10–15 cm 的 H1、O16；规则 11/11 出数前写下）：H1 z = %+.2f，O16 z = %+.2f → 【%s】" % (zs["H1"], zs["O16"], v))
    print("把上面整段输出发给 Claude。不管是哪一种结果，都照实写进 M3 小结后记和技术报告 6.4。")
