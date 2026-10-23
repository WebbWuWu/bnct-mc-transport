# -*- coding: utf-8 -*-
# ============================================================
#  check_kerma.py  --  M2：径迹长度 × MT=301 的 kerma 计数 + γ 产额计数，判据 K0–K4
#  判据归 Claude。2026-10-23 晚写（讲义第 11 课 §8 的 K1–K3）。被测：run_slab 新加的 6 个键（你写）。
#  运行：Windows，Ctrl+F5。水模体（build_water）+ 隐式俘获，跑三次：KUR 源 3000 条（K0、K2、K3）、
#        单能 1 MeV 2000 条（K1）、单能 1 MeV 1000 条「故意改坏」（K4）。约 1–2 分钟。全过最后打印「kerma 判据全过」。
# ============================================================
#
#  记号
#    m_k、SE_k   第 k 层的均值、标准误（都按「每源中子」，单位 eV 或「个」）
#    R_k         相对误差 SE_k / m_k。🔴 规则 14：R > 0.2 的层 SE 不可信，不拿去算 z
#    E_R         H(n,γ) 里氘核的反冲能（γ 带走 2.2233 MeV 的动量，氘核就要往反方向退一点）
#    f_801       B10 吸收里走 801 那一支（放 0.4778 MeV γ）的份额
#
#  判据（5 条，每条一行「抓什么」）
#    K0 结构：6 个新键都在；kerma 按核素名分、每条 NL 层；γ 两条各 NL 层；γ 个数处处 ≥ 0
#       抓：finish 没交出去、按核素 / 不按核素搞反（10/23 你写第 2 步时出过）、sq 被「换」成列表
#    K1 H1：径迹 kerma − E_R × gH（扣掉氘核反冲）  vs  碰撞估计器 sigma_x["H1"]（只有弹性反冲）
#       两条路：MT=301 查表 × 径迹（新）  vs  每次散射记 E − E′（第 8 课，10/17 和 OpenMC 判过）
#       比的是同一个量：H1 弹性反冲交给质子的能量（先扣掉 H1 MT=301 里多出来的氘核反冲 —— 纪律 2）
#       🔴 碰撞估计器很「吵」：KUR 源里快中子只占 6%，N = 4000 时连 5 层一块的 R 都 > 0.2（Claude 10/23 实测）。
#          所以 K1 单独用【单能 1 MeV】源跑 2000 条（每个中子都是快中子），并且每 5 层并成一块比：
#          块的 SE 取 5 层 SE 直接相加（层与层完全相关时的上限 —— 只会把 SE 估大、不会误报）
#          门槛：两边 R_块 ≤ 0.2 的块 |z| ≤ 3
#       分辨率有限（大约只抓 30% 以上的错）：抓的是接线错 —— heat 接到别的核素、少乘 BARN / 数密度、
#       fly 两处漏记一处、end_history 漏加。D_H 的精细比较靠 10/24–25 和 OpenMC 的 heating 逐层比
#    K2 B10：径迹 kerma  vs  逐道 Q 值的 Br（第 8 课）
#       两条路：MT=301（NJOY 算的「每次反应交给带电粒子多少」）  vs  σ800·Q800 + σ801·Q801
#       热区两者只差 B10 的弹性反冲和 600 / 700 两个小道；同一段路记的，涨落几乎全消掉
#       门槛（看数据前定）：R ≤ 0.2 的层 |kerma / Br − 1| ≤ 0.01 —— 我的预测是 1.000–1.002
#       抓：B10 的系数用错数密度 / 用了别的核素的 heat
#    K3 γ 产额：两条线各自核一次份额
#       ① gH / ag：ag 是所有核素 (n,γ) 之和，H1 的应占 ≥ 99.9%（O16 的 102 只有它的 0.03%，B10 的 102 约 0.002%）
#       ② gB × Q_val / Bc = 801 的份额 f_801（Bc = N_B (σ800+σ801) × Q_val，gB = N_B σ801；两个道都是 1/v，比值几乎不随能量变）
#          f_801 从库里另算一次（0.0253 eV 处 σ801 / (σ800 + σ801)），10/9 实测 0.93702；门槛 |差| ≤ 0.002
#       抓：gH 用成别的核素 / 别的道、gB 用成 800、少乘数密度
#    K4 试响（纪律 4）：把 H1 的 heat 乘 2 再跑一次（一个真的写错的输入：同一项记了两遍），K1 必须响；
#       另外把结果里的 gH 乘 1.01、B10 kerma 乘 1.02、gB 换成 800 的份额（假数据），K3①、K2、K3② 必须各自响
#       抓：门槛松到什么都放过
#
#  另外只打印、不判：O16 的径迹 kerma / 碰撞 D_H（差的就是「质心系各向同性」近似 #4 —— 10/25 方案乙的同一个问题）；
#                   γ 总能量（每源中子）
# ============================================================

import math
import os
import sys
import copy

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from lesson07_freegas import run_slab, Q_val
from lesson08_dose import kur_table, KEEP_FACTOR, D, NL
from lesson11_water import build_water

N_RUN = 3000
N_K1 = 2000
E_K1 = 1.0e6
SEED = 2310
N_FIRE = 1000
E_GAMMA_H = 2.2233e6       # eV，H(n,γ)²H 的 γ（库里 2223300 eV 那一条，讲义 §1.1）
E_GAMMA_B = 0.4778e6       # eV，⁷Li* → ⁷Li 的 γ
M_D_C2 = 1875.613e6        # eV，氘核静止能量（K1 的 E_R 第二条路）
Z_MAX = 3.0
R_MAX = 0.2
K2_TOL = 0.01
K3_H_MIN = 0.999
K3_B_TOL = 0.002
NEW_KEYS = ["kerma_sum", "kerma_sq", "gH_sum", "gH_sq", "gB_sum", "gB_sq"]


def mean_se(s, sq, N):
    """第 4 课：每层均值和标准误（每源中子）。"""
    ms, ses = [], []
    for k in range(len(s)):
        m = s[k] / N
        var = (sq[k] / N - m * m) * N / (N - 1)
        ms.append(m)
        ses.append(math.sqrt(max(var, 0.0) / N))
    return ms, ses


def run(nuclides, order, edges, shares, N, seed):
    return run_slab(N, seed, D, nuclides, None, order, z_src=0.0, implicit=True, NL=NL, src=(edges, shares))


def run_mono(nuclides, order, E0, N, seed):
    return run_slab(N, seed, D, nuclides, E0, order, z_src=0.0, implicit=True, NL=NL)


def deuteron_recoil(nuclides):
    """E_R 两条路：① 库里 H1 在表下限 E[0] 处 heat / σ₁₀₂（那里几乎全是俘获）② E_γ² / (2 M_d c²)。"""
    h = [n for n in nuclides if n["name"] == "H1"][0]
    i102 = h["mts"].index(102)
    lib = h["heat"][0] / h["XS"][i102][0]
    formula = E_GAMMA_H ** 2 / (2 * M_D_C2)
    return lib, formula


def f801_lib(nuclides):
    b = [n for n in nuclides if n["name"] == "B10"][0]
    i0 = b["mts"].index(800)
    i1 = b["mts"].index(801)
    from xs_lookup import sigma_at
    sig, _ = sigma_at(b["E"], b["XS"], 0.0253)
    return sig[i1] / (sig[i0] + sig[i1])


def k0(r, names):
    for key in NEW_KEYS:
        assert key in r, "K0：返回字典里没有 %r —— finish 没交出去？现有的键：%r" % (key, sorted(r))
    for key in ["kerma_sum", "kerma_sq"]:
        assert isinstance(r[key], dict) and sorted(r[key]) == sorted(names), \
            "K0：%s 应是「核素名 → 30 层」的字典，键应为 %r，实际是 %s" % (key, names, type(r[key]).__name__ if not isinstance(r[key], dict) else sorted(r[key]))
        for nm in names:
            assert len(r[key][nm]) == NL, "K0：%s[%r] 有 %d 层，应为 %d" % (key, nm, len(r[key][nm]), NL)
    for key in ["gH_sum", "gH_sq", "gB_sum", "gB_sq"]:
        assert isinstance(r[key], list) and len(r[key]) == NL, \
            "K0：%s 应是一条 %d 层的列表（γ 不按核素分，框架第 2 问 甲），实际是 %s" % (key, NL, type(r[key]).__name__)
        assert min(r[key]) >= 0.0, "K0：%s 有负数（个数不可能为负）：最小 %r" % (key, min(r[key]))
    assert max(r["gH_sum"]) > 0 and max(r["gB_sum"]) > 0 and max(max(r["kerma_sum"][nm]) for nm in names) > 0, \
        "K0：新计数全是 0 —— fly 里没记账，或者 end_history 没入账"


def k1(r, N, E_R, nucl="H1", block=5):
    """返回 (逐块的行, 没过的块号)。块 = 连续 block 层；块的 SE 取各层 SE 之和（完全相关的上限，偏保守）。"""
    kA, sA = mean_se(r["kerma_sum"][nucl], r["kerma_sq"][nucl], N)
    kC, sC = mean_se(r["sigma_x"][nucl], r["sigma_x2"][nucl], N)
    if nucl == "H1":
        g, sg = mean_se(r["gH_sum"], r["gH_sq"], N)
        kA = [kA[k] - E_R * g[k] for k in range(NL)]
        sA = [sA[k] + E_R * sg[k] for k in range(NL)]
    rows, bad = [], []
    for b0 in range(0, NL, block):
        ks = range(b0, min(b0 + block, NL))
        a = sum(kA[k] for k in ks); c = sum(kC[k] for k in ks)
        ea = sum(sA[k] for k in ks); ec = sum(sC[k] for k in ks)
        Ra = ea / abs(a) if a != 0 else float("inf")
        Rc = ec / abs(c) if c != 0 else float("inf")
        ok = Ra <= R_MAX and Rc <= R_MAX
        z = (a - c) / math.sqrt(ea ** 2 + ec ** 2) if ok else float("nan")
        rows.append((b0, a, c, a / c if c else float("nan"), z, Rc))
        if ok and abs(z) > Z_MAX:
            bad.append(b0 // block + 1)
    return rows, bad


def k2(r, N):
    kB, sB = mean_se(r["kerma_sum"]["B10"], r["kerma_sq"]["B10"], N)
    Br, sBr = mean_se(r["Br_sum"], r["Br_sq"], N)
    ratios, bad = [], []
    for k in range(NL):
        R = sBr[k] / Br[k] if Br[k] > 0 else float("inf")
        q = kB[k] / Br[k] if Br[k] > 0 else float("nan")
        ratios.append((k, q, R))
        if R <= R_MAX and abs(q - 1) > K2_TOL:
            bad.append(k + 1)
    return ratios, bad


def k3(r, N, f801):
    bad_H, bad_B, rows = [], [], []
    for k in range(NL):
        ag = r["ag_sum"][k]
        qH = r["gH_sum"][k] / ag if ag > 0 else float("nan")
        Bc = r["Bc_sum"][k]
        qB = r["gB_sum"][k] * Q_val / Bc if Bc > 0 else float("nan")
        rows.append((k, qH, qB))
        if ag > 0 and not (K3_H_MIN <= qH <= 1.0 + 1e-12):
            bad_H.append(k + 1)
        if Bc > 0 and abs(qB - f801) > K3_B_TOL:
            bad_B.append(k + 1)
    return rows, bad_H, bad_B


if __name__ == "__main__":
    print("=" * 76)
    edges, shares = kur_table()
    nuclides, order = build_water(KEEP_FACTOR * edges[-1])
    names = [n["name"] for n in nuclides]
    print("材料：%s；N = %d，种子 %d，隐式俘获，KUR 源，%d 层 × %.2f cm" % (names, N_RUN, SEED, NL, D / NL))

    E_R, E_R_f = deuteron_recoil(nuclides)
    print("E_R（氘核反冲）：库里 heat/σ₁₀₂ 在 E[0] = %.1f eV；公式 E_γ²/(2M_d c²) = %.1f eV" % (E_R, E_R_f))
    assert abs(E_R / E_R_f - 1) < 0.005, "E_R 两条路差超过 0.5%%：%.1f vs %.1f —— heat 接错了核素？" % (E_R, E_R_f)
    f801 = f801_lib(nuclides)
    print("f_801（801 份额，0.0253 eV）= %.5f（10/9 实测 0.93702）" % f801)

    r = run(nuclides, order, edges, shares, N_RUN, SEED)

    # ---------------- K0 ----------------
    k0(r, names)
    print("\nK0 过：6 个新键都在，结构对，γ 个数无负数")

    # ---------------- K1 ----------------
    r1 = run_mono(nuclides, order, E_K1, N_K1, SEED)
    rows, bad1 = k1(r1, N_K1, E_R)
    print("\nK1  H1：径迹 kerma − E_R·gH  vs  碰撞 D_H（单能 %.0f MeV × %d 条；eV / 源中子，每 5 层一块）" % (E_K1 / 1e6, N_K1))
    print("   块  深度cm         径迹（扣氘核）     碰撞        比值      z     碰撞R块")
    for b0, a, c, q, z, Rc in rows:
        print("  %3d  %5.1f–%-5.1f  %12.5g  %12.5g  %8.4f  %6.2f  %6.3f" % (b0 // 5 + 1, b0 * D / NL, (b0 + 5) * D / NL, a, c, q, z, Rc))
    nz = [z for _, _, _, _, z, _ in rows if not math.isnan(z)]
    assert nz, "K1：没有一块的 R ≤ 0.2 —— N 太小，告诉 Claude 加大 N_RUN"
    print("   参与判的块 %d / %d，max|z| = %.2f（门槛 %.0f）" % (len(nz), len(rows), max(abs(z) for z in nz), Z_MAX))
    assert not bad1, "K1 没过：第 %s 块 |z| > %.0f —— 先看这几块的比值是偏大还是偏小" % (bad1, Z_MAX)
    print("K1 过")

    # ---------------- K2 ----------------
    ratios, bad2 = k2(r, N_RUN)
    qs = [q for _, q, R in ratios if R <= R_MAX]
    print("\nK2  B10：径迹 kerma / 逐道 Q 值 Br（每层比值）：最小 %.5f  最大 %.5f（%d 层参与；预测 1.000–1.002）" % (min(qs), max(qs), len(qs)))
    print("    第 1、10、20、30 层：%s" % ", ".join("%.5f" % ratios[k][1] for k in (0, 9, 19, 29)))
    assert not bad2, "K2 没过：第 %s 层 |比值 − 1| > %.2f" % (bad2, K2_TOL)
    print("K2 过")

    # ---------------- K3 ----------------
    rows3, bH, bB = k3(r, N_RUN, f801)
    qHs = [q for _, q, _ in rows3 if not math.isnan(q)]
    qBs = [q for _, _, q in rows3 if not math.isnan(q)]
    print("\nK3① gH / ag：最小 %.6f  最大 %.6f（应在 [%.3f, 1]）" % (min(qHs), max(qHs), K3_H_MIN))
    print("K3② gB·Q_val / Bc：最小 %.5f  最大 %.5f（应在 f_801 = %.5f ± %.3f）" % (min(qBs), max(qBs), f801, K3_B_TOL))
    assert not bH, "K3① 没过：第 %s 层 gH 占 (n,γ) 总数的份额不在 [%.3f, 1] —— gH 用错核素 / 道 / 数密度？" % (bH, K3_H_MIN)
    assert not bB, "K3② 没过：第 %s 层 gB·Q_val/Bc 离 f_801 超过 %.3f —— gB 用成 800？" % (bB, K3_B_TOL)
    print("K3 过")

    # ---------------- 只打印 ----------------
    rowsO, _ = k1(r1, N_K1, E_R, nucl="O16")
    print("\n（只看不判）O16（单能 1 MeV 那次）：径迹 kerma / 碰撞 D_H —— 差的是近似 #4（质心系各向同性），10/25 方案乙同一个问题")
    print("    各块（每 5 层）：%s" % ", ".join("%.3f（R=%.2f）" % (row[3], row[5]) for row in rowsO))
    gHm, _ = mean_se(r["gH_sum"], r["gH_sq"], N_RUN)
    gBm, _ = mean_se(r["gB_sum"], r["gB_sq"], N_RUN)
    EH = sum(gHm) * E_GAMMA_H
    EB = sum(gBm) * E_GAMMA_B
    print("（只看不判）每个源中子在 15 cm 里产生的 γ 能量：H 线 %.4g MeV + B 线 %.4g MeV = %.4g MeV"
          % (EH / 1e6, EB / 1e6, (EH + EB) / 1e6))

    # ---------------- K4 试响 ----------------
    print("\nK4 试响（纪律 4）")
    bad_nuc = copy.deepcopy(nuclides)
    iH = names.index("H1")
    bad_nuc[iH]["heat"] = nuclides[iH]["heat"] * 2.0        # 一个真的写错的输入：H1 的 heat 记了两遍（截面不动，随机数序列不变）
    rb = run_mono(bad_nuc, order, E_K1, N_FIRE, SEED)
    _, b1 = k1(rb, N_FIRE, E_R)
    assert b1, "K4：H1 的 heat × 2，K1 居然没响 —— 门槛太松"
    print("  H1 heat × 2 → K1 在 %d 块响 ✓" % len(b1))
    r2 = copy.deepcopy(r)
    r2["gH_sum"] = [x * 1.01 for x in r2["gH_sum"]]
    _, bH2, _ = k3(r2, N_RUN, f801)
    assert bH2, "K4：gH × 1.01，K3① 没响"
    print("  gH × 1.01 → K3① 在 %d 层响 ✓" % len(bH2))
    r3 = copy.deepcopy(r)
    r3["kerma_sum"]["B10"] = [x * 1.02 for x in r3["kerma_sum"]["B10"]]
    _, b2 = k2(r3, N_RUN)
    assert b2, "K4：B10 kerma × 1.02，K2 没响"
    print("  B10 kerma × 1.02 → K2 在 %d 层响 ✓" % len(b2))
    r4 = copy.deepcopy(r)
    r4["gB_sum"] = [x * (1 - f801) / f801 for x in r4["gB_sum"]]      # 相当于数成了 800 那一支
    _, _, bB4 = k3(r4, N_RUN, f801)
    assert bB4, "K4：gB 换成 800 的份额，K3② 没响"
    print("  gB 换成 800 那一支 → K3② 在 %d 层响 ✓" % len(bB4))

    print("\nkerma 判据全过（K0–K4）")
