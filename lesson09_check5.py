# -*- coding: utf-8 -*-
# ============================================================
#  lesson09_check5.py  --  第 9 课判据 5：OpenMC 加 / 不加 S(α,β)，量出「自由气体近似」差多少
#  判据归 Claude，2026-10-16 21:20 写。被测：outputs\ 下 4 份 OpenMC 结果（都是 100 批 × 1 万）
#     自由气体：openmc_kur.json（种子 1）、openmc_kur_s2.json（种子 2）
#     加 S(α,β)：openmc_kur_sab3.json（种子 3）、openmc_kur_sab4.json（种子 4）—— 只给 H1 加了 c_H_in_H2O
#  运行：Windows，Ctrl+F5（只用 math / json / os）。几秒。
# ============================================================
#
#  这一条在干什么（讲义 §9 判据 5）
#    不是「对不对」，是「近似差多少」：你的程序用的是自由气体（近似 #9），真实组织里的氢是被化学键
#    绑住的，正确做法是 S(α,β) 热散射律。两边都用 OpenMC 算，只差这一张表，差出来的就是近似 #9 的大小。
#    —— 这是面试材料：「我知道我的程序哪里不对、不对多少」。
#
#  记号
#    「合并」   同一设置两次独立运行取平均：m = (m₁ + m₂)/2，SE = sqrt(SE₁² + SE₂²)/2（和判据 3 一样）
#    相对差     (S(α,β) − 自由气体) / 自由气体，按层
#    z          (S(α,β) − 自由气体) / sqrt(SE² + SE²)，按层（四次运行种子都不同，互相独立）
#    峰深       每层的值取最大的那层，再用它和左右两层拟一条抛物线，取顶点的深度（比 0.5 cm 一格细）
#
#  判据（每条旁边写抓什么）
#    5-1 S(α,β) 真的加上了、也只加在该加的那两次：
#        kur_sab3、kur_sab4 的 materials.xml 里有 <sab name="c_H_in_H2O"/>；kur、kur_s2 里没有
#        抓：表名拼错（H20 ← 数字 0）OpenMC 会报错，但「加在了循环外面 / 加到了自由气体那次 / 第二圈忘了加」
#            都不报错 —— 只能看 OpenMC 实际读到的输入文件
#    5-2 两次 S(α,β) 互比（种子 3 vs 4），和判据 3-5 同样的门槛：χ² ∈ [13.79, 53.67]、max|z| ≤ 4
#        只判 4 个主量（通量、D_B、D_N、(n,γ)）；3 个 heating 只打印
#        抓：两次其实不是同一个设置（比如一次加了表一次没加、两圈写进了同一个目录又被覆盖）
#        ⚠️ 10/16 实测：O16 heating 两次互比 χ² = 74（门槛外），H1 / C12 heating 48 / 51（边缘）；
#           自由气体那两次 O16 heating 是 26。heating 由少数快中子碰撞主导，批均值分布拖尾，
#           OpenMC 的批间 SE 可能估小 —— 和 S(α,β) 无关（heating 前后只差 < 0.4%），但周六判据 4（D_H）要先查它。
#    5-3 只量、不判：4 个量（+ 3 个 heating）逐层相对差和 z、30 层总和之比、D_B 峰深前后对比
#        —— 对你 10/16 17:25 的猜测：「D_B 峰往浅处移（S(α,β) 让热区散射截面变大，热中子扩散得近）」
# ============================================================

import math
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "outputs")
RUN_DIR = os.path.join(OUT_DIR, "openmc")

FREE = [("自由气体 种子1", "openmc_kur.json", "kur"), ("自由气体 种子2", "openmc_kur_s2.json", "kur_s2")]
SAB = [("S(a,b) 种子3", "openmc_kur_sab3.json", "kur_sab3"), ("S(a,b) 种子4", "openmc_kur_sab4.json", "kur_sab4")]
SAB_TAG = '<sab name="c_H_in_H2O"'

NL = 30
DZ = 0.5
Q_B = 2341569.0          # eV，同判据 3（只用来把 D_B 换成 eV 打印）

NAMES = [("flux", "kur flux"), ("D_B", "kur B10 n,a"), ("D_N", "kur N14 n,p"), ("(n,g)", "kur n,gamma"),
         ("H1 heat", "kur H1 heating"), ("C12 heat", "kur C12 heating"), ("O16 heat", "kur O16 heating")]
CHI2_LO, CHI2_HI = 13.79, 53.67
ZMAX = 4.0


def load(fname):
    p = os.path.join(OUT_DIR, fname)
    assert os.path.exists(p), "找不到 %s —— 先在 WSL 里跑 lesson09_openmc.py（run_Sab = True）" % p
    with open(p, "r") as f:
        return json.load(f)


def merge(a, b, name):
    m = [(x + y) / 2 for x, y in zip(a[name]["mean"], b[name]["mean"])]
    s = [math.sqrt(x * x + y * y) / 2 for x, y in zip(a[name]["se"], b[name]["se"])]
    return m, s


def peak_depth(m):
    """最大的那层 + 左右两层拟抛物线，返回 (第几层（从 1 数）, 顶点深度 cm)。"""
    k = max(range(NL), key=lambda i: m[i])
    if k == 0 or k == NL - 1:
        return k + 1, (k + 0.5) * DZ
    y0, y1, y2 = m[k - 1], m[k], m[k + 1]
    off = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)       # 以层为单位，相对第 k 层中心
    return k + 1, (k + 0.5 + off) * DZ


if __name__ == "__main__":
    # ---------------- 5-1 ----------------
    print("=" * 78)
    print("5-1 S(α,β) 加没加：读 OpenMC 实际用的 materials.xml")
    for lab, _, d in FREE + SAB:
        p = os.path.join(RUN_DIR, d, "materials.xml")
        assert os.path.exists(p), "找不到 %s" % p
        with open(p, "r") as f:
            has = SAB_TAG in f.read()
        want = (d.startswith("kur_sab"))
        print("  %-16s %-10s 有 c_H_in_H2O：%s" % (lab, d, "是" if has else "否"))
        assert has == want, "%s（%s）：materials.xml 里%s c_H_in_H2O，应该%s" % (lab, d, "有" if has else "没有", "有" if want else "没有")
    print("  5-1 过")

    f1, f2 = load(FREE[0][1]), load(FREE[1][1])
    s3, s4 = load(SAB[0][1]), load(SAB[1][1])
    for o in (f1, f2, s3, s4):
        assert o.get("batches") == 100 and o.get("particles") == 10000, "有一份不是 100 批 × 1 万"

    # ---------------- 5-2 ----------------
    print("=" * 78)
    print("5-2 两次 S(α,β) 互比（种子 3 vs 4）：χ² 应 ≈ 30")
    for disp, name in NAMES:
        judged = not name.endswith("heating")
        zs = [(s3[name]["mean"][k] - s4[name]["mean"][k]) /
              math.sqrt(s3[name]["se"][k] ** 2 + s4[name]["se"][k] ** 2) for k in range(NL)]
        chi2 = sum(z * z for z in zs)
        zmax = max(abs(z) for z in zs)
        print("  %-10s χ² %6.1f   max|z| %5.2f   %s" % (disp, chi2, zmax, "" if judged else "（只看，不判；见文件头 ⚠️）"))
        assert (not judged) or CHI2_LO <= chi2 <= CHI2_HI and zmax <= ZMAX, \
            "%s：两次 S(α,β) 互相对不上（χ² %.1f、max|z| %.2f）—— 两次不是同一个设置？" % (disp, chi2, zmax)
    print("  5-2 过（4 个主量）")

    # ---------------- 5-3 ----------------
    print("=" * 78)
    print("5-3 量：S(α,β) 合并 vs 自由气体合并（只量、不判）")
    print("  %-10s %10s   相对差（%%），第 1、4、7、…、28 层" % ("量", "Σ30层之比"))
    for disp, name in NAMES:
        F = merge(f1, f2, name)
        S = merge(s3, s4, name)
        rel = [S[0][k] / F[0][k] - 1 for k in range(NL)]
        print("  %-10s %10.4f   %s" % (disp, sum(S[0]) / sum(F[0]), " ".join("%+5.1f" % (100 * r) for r in rel[::3])))
    print()
    F = merge(f1, f2, "kur B10 n,a")
    S = merge(s3, s4, "kur B10 n,a")
    print("  D_B 逐层（eV / 源中子；15 ppm）")
    print("  %4s %8s %12s %12s %8s %7s" % ("层", "深 cm", "自由气体", "S(α,β)", "相对差", "z"))
    for k in range(NL):
        z = (S[0][k] - F[0][k]) / math.sqrt(S[1][k] ** 2 + F[1][k] ** 2)
        print("  %4d %8.2f %12.1f %12.1f %+7.1f%% %+7.1f" % (k + 1, (k + 0.5) * DZ, F[0][k] * Q_B, S[0][k] * Q_B,
                                                             100 * (S[0][k] / F[0][k] - 1), z))
    kf, df = peak_depth(F[0])
    ks, ds = peak_depth(S[0])
    print()
    print("  D_B 峰：自由气体 第 %d 层，抛物线顶点 %.2f cm；S(α,β) 第 %d 层，顶点 %.2f cm；移动 %+.2f cm"
          % (kf, df, ks, ds, ds - df))
    print("  D_B 峰值：自由气体 %.1f eV → S(α,β) %.1f eV（%+.1f%%）"
          % (max(F[0]) * Q_B, max(S[0]) * Q_B, 100 * (max(S[0]) / max(F[0]) - 1)))
    print("=" * 78)
    print("判据 5：5-1、5-2 过；5-3 的数进偏差表「已知近似 #9」一栏")
