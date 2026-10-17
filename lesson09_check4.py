# -*- coding: utf-8 -*-
# ============================================================
#  lesson09_check4.py  --  第 9 课判据 4：D_H（弹性反冲）按核素和 OpenMC 的 heating 比
#  判据归 Claude，2026-10-17 写。被测：你的 run_slab 新加的「按核素分开的 D_H」（框架第一轮第 5 问甲 + 10/17 第 1 问乙）
#  运行：Windows，Ctrl+F5（4-1 要 h5py 读库；其余只用 math / json / os）。几秒。
#  🔴 跑之前改文件头下面「要你填的三个常数」—— 等第 1、2、3 问答完我替你填好
# ============================================================
#
#  讲义 §7：两边不是同一个量，要先让它们变成同一个量
#    你的 D_H     每次弹性散射 wt × (E − E')，按撞到的核素记（自由气体；热区上散射为负）
#    OpenMC       score 'heating' 按核素 = 库里的 KERMA（MT=301）× 径迹长度。KERMA 包含该核素【所有】反应留下的能量
#    方案甲：只比 H1、C12、O16（源 ≤ 2 MeV 时，它们除了弹性只有 (n,γ)）
#      H1：KERMA 里还有 H(n,γ) 的氘核反冲 —— 每次约 1.3 keV，每层占 H1 heating 的 0.6–3.3%（10/17 实测），
#          远大于 SE，必须扣：OpenMC H1 heating − E_R × (n,γ)。E_R 由 4-1 从库里量出来
#      C12、O16：(n,γ) 次数只有 H 的千分之一量级（σ_γ 0.0035 / 0.0002 b vs H 0.332 b），反冲不扣，偏差 < 0.2%
#      ⚠️ OpenMC 的 '(n,gamma)' tally 是 5 个核素加起来的（H 占 ~99.4%），按全部算 H 的 → 多扣 ~0.6% × 3% ≈ 0.02%，忽略
#    §7 ②「KERMA 假设靶核静止」：热区弹性那部分 KERMA ≈ 0.26 eV·b，H 俘获那部分 ≈ 438 eV·b，前者占 0.06% → 忽略
#
#  记号
#    合并（你）   两个种子：m = (m₁+m₂)/2，SE = sqrt(SE₁²+SE₂²)/2（判据 3 同）
#    合并（OMC） OMC_RUNS 里 M 次：m = 平均，SE = sqrt(ΣSE²)/M
#    z_k        (你的 − OpenMC) / sqrt(SE你² + SE_OMC²)，第 k 层
#    E_R        H(n,γ) 一次留下的反冲能（eV）
#
#  判据（≤ 5 条，每条旁边写抓什么）
#    4-1 库里的 KERMA 扣法对不对：两条独立的路算 E_R
#        路 1（库）：热区（0.005–1 eV）每个格点 (KERMA − σ_el·E·2A/(A+1)²) / σ_γ，取中位数
#        路 2（公式）：Q² / (2·(A+1)·m_n c²)，Q = 库里 MT=102 的 Q 值
#        门槛：差 < 3%
#        抓：KERMA 里其实含了 γ 能量（会大 ~1700 倍）或根本没算俘获反冲（≈ 0）—— 那样扣法就是错的
#    4-2 回归：只证明「没改坏」（纪律 4）
#        (a) 新结果的每一个旧键，和改之前同种子的旧 json 逐位相同（两个种子都比）
#            抓：改动碰到了随机数顺序或别的打分（缩进挪了一行、多抽了一次 rng.random()……）
#        (b) 新字典 5 个核素逐层加起来 == 旧的 DH_sum（相对 1e-9）；键正好是 5 个核素名
#            抓：漏了某个核素 / 记错层（k0 vs k_hit）/ 只在一个分支里记了 / 历史结束忘了累加
#    4-3 真比：H1（扣反冲后）、C12、O16 各 30 层（门槛按第 3 问，见 MODE）
#        抓：两边物理不同（含讲义第 6 课近似 4「CM 各向同性」—— 10/17 粗看深处你偏低 15–25%）
#    4-4 判据自己会响：(a) 把 H1 某一层挪到隔壁层 → 4-2(b) 必须响；(b) H1 和 O16 标签对调 → 4-3 必须响
#        抓：判据太松、什么都放过
#    4-5 只量、不判：按深度 5 段（每段 3 cm）列「你 / OpenMC」比值 ± 误差，三个核素 + 反冲占比
#        —— 给偏差表 D_H 那一行写「已知原因」用
# ============================================================

import math
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "outputs")

# ==================== 要你填的三个常数（第 1、2、3 问答完由我填）====================
KEY_SUM = "sigma_x"           # ⑤ 你起的键名（10/17 用户写）：按核素的 Σx（字典：核素名 → 30 层）
KEY_SQ = "sigma_x2"           # ⑤ 你起的键名（10/17 用户写）：按核素的 Σx²
NEW_STEM = "dose_N100000_s%d_nuc"     # 第 2 问甲（10/17 定）：重跑后的新文件名（不带 .json），%d = 种子
MODE = "甲"                   # 第 3 问：判据 4-3 的门槛（见下）
# =================================================================================

OLD_STEM = "dose_N100000_s%d"
SEEDS = [1010, 1200]
N_HIST = 100000
NL = 30
DZ = 0.5
NUC5 = ["H1", "C12", "N14", "O16", "B10"]
CMP = [("H1", "kur H1 heating"), ("C12", "kur C12 heating"), ("O16", "kur O16 heating")]
OMC_RUNS = ["openmc_kur.json", "openmc_kur_s2.json"] + ["openmc_kur_s%d.json" % s for s in range(5, 13)]
M_N_EV = 939.56542052e6       # 中子静能 eV（CODATA 2018）
ER_TOL = 0.03
ZMAX = 4.0
CHI2_LO, CHI2_HI = 13.79, 53.67   # 30 自由度 99%（只在 MODE 丙 用）
#  MODE 的三种（第 3 问）：
#    甲  max|z| ≤ 4，且 χ² ≤ 「OpenMC 自己 45 对互比里最大的那个 χ²」（门槛从数据里标，已含层间相关）
#    乙  只判 max|z| ≤ 4（单层 z 不受层间相关影响）；χ² 只打印
#    丙  照搬判据 3：χ² ∈ [13.79, 53.67] 且 max|z| ≤ 4（check_se 已证明 heating 会误报 ~22%，不推荐）


def load(fname):
    p = os.path.join(OUT_DIR, fname)
    assert os.path.exists(p), "找不到 %s" % p
    with open(p, "r") as f:
        return json.load(f)


def my_mean_se(s, sq):
    """Σx、Σx² → 30 层均值、SE（第 4 课公式，和 mcstat.mean_se 同一个）。"""
    m, e = [], []
    for k in range(NL):
        mu = s[k] / N_HIST
        var = (sq[k] / N_HIST - mu * mu) * N_HIST / (N_HIST - 1)
        m.append(mu)
        e.append(math.sqrt(max(var, 0.0) / N_HIST))
    return m, e


def merge2(a, b):
    return ([(x + y) / 2 for x, y in zip(a[0], b[0])],
            [math.sqrt(x * x + y * y) / 2 for x, y in zip(a[1], b[1])])


def merge_omc(runs, name):
    M = len(runs)
    m = [sum(o[name]["mean"][k] for o in runs) / M for k in range(NL)]
    e = [math.sqrt(sum(o[name]["se"][k] ** 2 for o in runs)) / M for k in range(NL)]
    return m, e


def zstats(a, b):
    """a、b = (均值列表, SE 列表)。返回 z 列表、1σ 内、2σ 内、max|z|、χ²。"""
    z = [(a[0][k] - b[0][k]) / math.sqrt(a[1][k] ** 2 + b[1][k] ** 2) for k in range(NL)]
    return (z, sum(1 for x in z if abs(x) < 1), sum(1 for x in z if abs(x) < 2),
            max(abs(x) for x in z), sum(x * x for x in z))


def judge(st, chi2_cap):
    z, n1, n2, zmax, chi2 = st
    if MODE == "甲":
        return zmax <= ZMAX and chi2 <= chi2_cap
    if MODE == "乙":
        return zmax <= ZMAX
    return zmax <= ZMAX and CHI2_LO <= chi2 <= CHI2_HI


def er_from_library():
    """4-1：返回 (路 1 E_R 中位数, 路 1 的最小、最大, 路 2 E_R)。要 h5py。"""
    import h5py
    from xs_lookup import LIB
    p = os.path.join(LIB, "H1.h5")
    assert os.path.exists(p), "找不到 %s" % p
    with h5py.File(p, "r") as f:
        root = f["H1"]
        E = root["energy"]["294K"][:]
        A = float(root.attrs["atomic_weight_ratio"])

        def full(mt):
            ds = root["reactions"]["reaction_%03d" % mt]["294K"]["xs"]
            thr = int(ds.attrs["threshold_idx"])
            x = [0.0] * len(E)
            v = ds[:]
            for i in range(len(v)):
                x[thr + i] = float(v[i])
            return x

        el, cap, ker = full(2), full(102), full(301)
        Q = float(root["reactions"]["reaction_102"].attrs["Q_value"])
    fel = 2 * A / (A + 1) ** 2
    vals = [(ker[i] - el[i] * E[i] * fel) / cap[i] for i in range(len(E)) if 0.005 <= E[i] <= 1.0 and cap[i] > 0]
    assert len(vals) >= 5, "热区只找到 %d 个格点，读库的方式不对？" % len(vals)
    vals.sort()
    er_lib = vals[len(vals) // 2]
    er_formula = Q * Q / (2 * (A + 1) * M_N_EV)
    return er_lib, vals[0], vals[-1], er_formula, len(vals)


if __name__ == "__main__":
    # ---------------- 4-1 ----------------
    print("=" * 78)
    print("4-1 H(n,γ) 反冲能 E_R：库里的 KERMA vs 公式")
    er_lib, er_lo, er_hi, er_f, nv = er_from_library()
    print("  路 1（库，热区 %d 个格点）中位数 %.1f eV（范围 %.1f – %.1f）" % (nv, er_lib, er_lo, er_hi))
    print("  路 2（公式 Q²/(2(A+1)m_n c²)）    %.1f eV" % er_f)
    assert abs(er_lib / er_f - 1) < ER_TOL, \
        "库里 KERMA 量出的俘获反冲 %.1f eV 和公式 %.1f eV 差 %.1f%%（> 3%%）—— 扣法要改，先停" \
        % (er_lib, er_f, 100 * (er_lib / er_f - 1))
    E_R = er_lib
    print("  4-1 过：H1 heating 扣 %.1f eV × (n,γ)" % E_R)

    # ---------------- 4-2 ----------------
    print("=" * 78)
    print("4-2 回归：只证明没改坏")
    new = {}
    for sd in SEEDS:
        old = load(OLD_STEM % sd + ".json")
        nw = load(NEW_STEM % sd + ".json")
        new[sd] = nw
        diff = [k for k in old if old[k] != nw.get(k)]
        print("  种子 %d：旧键 %d 个，和旧 json 不同的 %d 个 %s" % (sd, len(old), len(diff), diff[:5]))
        assert not diff, "种子 %d：改完以后这些旧结果变了：%s —— 改动碰到了别的东西" % (sd, diff)
        assert KEY_SUM in nw and KEY_SQ in nw, "种子 %d：新 json 里没有键 %s / %s（⑤ 的键名对吗？）" % (sd, KEY_SUM, KEY_SQ)
        assert sorted(nw[KEY_SUM].keys()) == sorted(NUC5), \
            "种子 %d：按核素字典的键是 %s，应该正好是 %s" % (sd, sorted(nw[KEY_SUM].keys()), sorted(NUC5))
        for k in range(NL):
            tot = sum(nw[KEY_SUM][n][k] for n in NUC5)
            assert math.isclose(tot, nw["DH_sum"][k], rel_tol=1e-9, abs_tol=1e-9), \
                "种子 %d 第 %d 层：5 个核素加起来 %.9g ≠ 旧的 DH_sum %.9g" % (sd, k + 1, tot, nw["DH_sum"][k])
    print("  4-2 过：旧结果逐位不变；5 个核素逐层加起来 = DH_sum")

    # ---------------- 准备两边 ----------------
    omc = [load(f) for f in OMC_RUNS]
    ng = merge_omc(omc, "kur n,gamma")
    om = {}
    for nuc, tname in CMP:
        m, e = merge_omc(omc, tname)
        if nuc == "H1":
            m = [m[k] - E_R * ng[0][k] for k in range(NL)]
            e = [math.sqrt(e[k] ** 2 + (E_R * ng[1][k]) ** 2) for k in range(NL)]
        om[nuc] = (m, e)
    mine = {}
    for nuc, _ in CMP:
        a = my_mean_se(new[SEEDS[0]][KEY_SUM][nuc], new[SEEDS[0]][KEY_SQ][nuc])
        b = my_mean_se(new[SEEDS[1]][KEY_SUM][nuc], new[SEEDS[1]][KEY_SQ][nuc])
        mine[nuc] = merge2(a, b)
    caps = {}
    for nuc, tname in CMP:
        caps[nuc] = max(sum((omc[i][tname]["mean"][k] - omc[j][tname]["mean"][k]) ** 2 /
                            (omc[i][tname]["se"][k] ** 2 + omc[j][tname]["se"][k] ** 2) for k in range(NL))
                        for i in range(len(omc)) for j in range(i + 1, len(omc)))

    # ---------------- 4-3 ----------------
    print("=" * 78)
    print("4-3 真比：你的（%d 个种子合并）vs OpenMC（%d 次合并），MODE = %s" % (len(SEEDS), len(omc), MODE))
    print("  %-5s %6s %6s %7s %7s %12s %10s  判" % ("核素", "1σ内", "2σ内", "max|z|", "χ²", "χ² 上限(甲)", "Σ你/ΣOMC"))
    res = {}
    for nuc, _ in CMP:
        st = zstats(mine[nuc], om[nuc])
        ok = judge(st, caps[nuc])
        res[nuc] = ok
        print("  %-5s %6d %6d %7.2f %7.1f %12.1f %10.4f  %s" % (nuc, st[1], st[2], st[3], st[4], caps[nuc],
                                                                sum(mine[nuc][0]) / sum(om[nuc][0]), "过" if ok else "✗"))
        print("        z：" + " ".join("%+.1f" % x for x in st[0]))
    bad = [n for n in res if not res[n]]
    if bad:
        print("  ⚠️ 4-3 没过：%s —— 先看 4-5 的分段表，再定原因（不改门槛）" % bad)
    else:
        print("  4-3 过")

    # ---------------- 4-4 ----------------
    print("=" * 78)
    print("4-4 判据自己会响")
    nw = new[SEEDS[0]]
    shifted = list(nw[KEY_SUM]["H1"])
    shifted[5], shifted[6] = 0.0, shifted[6] + shifted[5]          # 第 6 层的 H1 记到第 7 层去
    rang = any(not math.isclose(sum((shifted if n == "H1" else nw[KEY_SUM][n])[k] for n in NUC5),
                                nw["DH_sum"][k], rel_tol=1e-9, abs_tol=1e-9) for k in range(NL))
    print("  (a) H1 第 6 层挪到第 7 层：4-2(b) %s" % ("响" if rang else "没响"))
    assert rang, "4-2(b) 抓不到记错层"
    st_sw = zstats(mine["O16"], om["H1"])
    rang2 = not judge(st_sw, caps["H1"])
    print("  (b) 你的 O16 当成 H1 去比：max|z| %.1f → 4-3 %s" % (st_sw[3], "响" if rang2 else "没响"))
    assert rang2, "4-3 抓不到核素标签对调"
    for f in (1.03, 0.97):
        st_f = zstats(([x * f for x in mine["H1"][0]], mine["H1"][1]), om["H1"])
        print("  （只看）你的 H1 整体 ×%.2f：max|z| %.2f、χ² %.1f → %s" % (f, st_f[3], st_f[4], "会响" if not judge(st_f, caps["H1"]) else "不响"))
    print("  4-4 过")

    # ---------------- 4-5 ----------------
    print("=" * 78)
    print("4-5 只量：按深度分 5 段的「你 / OpenMC」（± 是 1σ）")
    segs = [(0, 6), (6, 12), (12, 18), (18, 24), (24, 30)]
    print("  %-5s " % "核素" + " ".join("%17s" % ("%.0f–%.0f cm" % (a * DZ, b * DZ)) for a, b in segs))
    for nuc, _ in CMP:
        row = []
        for a, b in segs:
            sm = sum(mine[nuc][0][a:b]); so = sum(om[nuc][0][a:b])
            em = math.sqrt(sum(x * x for x in mine[nuc][1][a:b]))      # 段内层独立的估法，偏小（层间正相关）
            eo = math.sqrt(sum(x * x for x in om[nuc][1][a:b]))
            r = sm / so
            row.append("%7.3f ± %-7.3f" % (r, r * math.sqrt((em / sm) ** 2 + (eo / so) ** 2)))
        print("  %-5s " % nuc + " ".join("%17s" % s for s in row))
    rec = [E_R * ng[0][k] for k in range(NL)]
    h1raw = merge_omc(omc, "kur H1 heating")[0]
    print("  H1 反冲占 OpenMC H1 heating：" + " ".join("%.1f%%" % (100 * rec[k] / h1raw[k]) for k in range(0, NL, 5)) + "（第 1、6、…、26 层）")
    print("  ⚠️ ± 按段内层独立算，层间正相关时偏小，只看量级")
    print("=" * 78)
    print("判据 4：4-1、4-2、4-4 必须过；4-3 按 MODE %s；4-5 进偏差表 D_H 一行" % MODE)
