# -*- coding: utf-8 -*-
# ============================================================
#  check_kn_te.py  --  M4 判据 KT：康普顿「平均交出比例」⟨T⟩/E 对不对、做成表以后还对不对
#  判据归 Claude。2026-11-24 写。运行：Windows Ctrl+F5（openmc (3.10.20) 环境），要读光子库，几秒。
#  ⚠️ 写的时候 Claude 那边没有 h5py 和真库：KT-1、自检在复刻的函数上跑过；KT-2 用复刻的表跑过。
#     第一次在你电脑上跑的完整输出发给 Claude；它自己报错（不是 assert 响）多半是 Claude 的锅。
#
#  被测：photon.kn_T_E(E)          现算：μ 上 2001 点梯形积分，1 − ∫f·k dμ / ∫f dμ
#        photon.kn_TE_tab(pho, E)  查表：load_photon 里 1 keV–2.3 MeV 200 点等比，np.interp 线性插
#  先确认比的是同一个量（纪律 2）：
#    f(μ) = k²(k + 1/k − 1 + μ²) ∝ dσ_KN/dμ，k = E′/E。一次康普顿交给电子 T = E − E′，所以
#    ⟨T⟩/E = 1 − ⟨k⟩ = ∫(1 − k) dσ / ∫dσ = σ_tr / σ_KN    —— 左边是你算的，右边是教科书闭式解
#
#  三条：
#    [KT-1] 现算 vs 闭式解（σ_tr / σ_KN，Attix 式 7.13 / 7.11 那套；不经过 f(μ) 的数值积分）
#           6 个能量（光子会有的能量 + 表的两头附近，不取整数特殊值），相对差 < 10⁻³
#           抓：f(μ) 写错、少了「1 −」、积分分子漏乘 k、mₑc² 单位错
#    [KT-2] 表 vs 现算：1 keV–2.3 MeV 之间对数均匀随机 2000 个能量 + 表两头 + 两条 γ 线，最大相对差 < 10⁻³
#           抓：表的网格太粗、表错位（TE[i] 存的不是 E[i] 的值）、插值用错数组
#    [KT-3] 覆盖与越界：截断能量 E_CUT_G 和最高的 γ（2.2233 MeV）都在表里；表外的能量 assert 必须响
#           抓：表的范围比程序会用到的窄（以后加 14 MeV 源时这条会先响）；越界 assert 被删
#  门槛：10⁻³ 是 11/23 定的（写在提示词第六节），看数之前。
#    依据：2001 点梯形的误差 ~ h² ~ 10⁻⁶；200 点等比表线性插 ~ 10⁻⁴（11/23 Claude 粗算）；留一个数量级。
#  抓不到什么（照实写）：
#    · 纯 KN 本身对不对（康普顿用自由电子，近似 #32）—— 归 OpenMC 光子对拍
#    · μ_tr 里光电、电子对两项 —— 归 μ_tr vs NIST（photon.py __main__ 末尾）
#  自检（纪律 4）：6 个故意写错的版本必须响，正确的模板必须不响（证明模板本身没写错）
# ============================================================
import os, sys
_MY_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # my_code 根目录
for _sub in ("whiteboard", "lessons", "openmc_compare", "checks", ""):
    _p = os.path.join(_MY_CODE, _sub) if _sub else _MY_CODE
    if _p not in sys.path:
        sys.path.insert(0, _p)

import math
import time
import numpy as np

import photon as PH

GATE = 1.0e-3
ME_C2_REF = 0.51099895e6            # eV，CODATA 2018；判据这边自己再写一遍，不用 PH.ME_C2（那样 mₑc² 写错两边一起错）
E_KT1 = [1.37e3, 3.7e4, 1.18e5, PH.E_G_B, 1.31e6, PH.E_G_H]   # 1.37 keV · 37 keV · 118 keV · B10 的 γ · 1.31 MeV · H 俘获 γ
N_KT2, SEED_KT2 = 2000, 61403


def te_closed(E):
    """⟨T⟩/E 的闭式解 = σ_tr / σ_KN（两个都除掉了公共因子 2π rₑ²）。"""
    a = E / ME_C2_REF
    b = 1.0 + 2.0 * a
    L = math.log1p(2.0 * a)                                  # ln(1 + 2α)；低能时比 log(1+2α) 准
    s_kn = (1 + a) / a**2 * (2 * (1 + a) / b - L / a) + L / (2 * a) - (1 + 3 * a) / b**2
    s_tr = (2 * (1 + a)**2 / (a**2 * b) - (1 + 3 * a) / b**2
            - (1 + a) * (2 * a * a - 2 * a - 1) / (a**2 * b**2) - 4 * a * a / (3 * b**3)
            - ((1 + a) / a**3 - 1 / (2 * a) + 1 / (2 * a**3)) * L)
    return s_tr / s_kn


# ---------------- 三条判据 ----------------
def kt1(fn=PH.kn_T_E, verbose=True):
    ok = True
    for E in E_KT1:
        got, ref = fn(E), te_closed(E)
        rel = got / ref - 1.0 if ref != 0 else float("inf")
        good = abs(rel) < GATE
        ok = ok and good
        if verbose:
            print("    E = %-10.5g eV  现算 %.6f  闭式解 %.6f  相对差 %+.2e  %s" % (E, got, ref, rel, "✅" if good else "❌"))
    return ok


def kt2_energies(pho):
    Eg = pho["KN"]["E"]
    rng = np.random.default_rng(SEED_KT2)
    Es = np.exp(rng.uniform(math.log(Eg[0]), math.log(Eg[-1]), N_KT2))
    return [float(Eg[0]), float(Eg[-1]), PH.E_G_B, PH.E_G_H] + [float(x) for x in Es]


def kt2(pho, tab=PH.kn_TE_tab, exact=PH.kn_T_E, verbose=True):
    worst, E_worst = 0.0, None
    for E in kt2_energies(pho):
        rel = abs(tab(pho, E) / exact(E) - 1.0)
        if rel > worst:
            worst, E_worst = rel, E
    ok = worst < GATE
    if verbose:
        print("    %d 个能量，最大相对差 %.2e（在 E = %.5g eV）  门槛 %.0e  %s"
              % (N_KT2 + 4, worst, E_worst, GATE, "✅" if ok else "❌"))
    return ok


def kt3(pho, verbose=True):
    Eg = pho["KN"]["E"]
    lines, ok = [], True
    c1 = Eg[0] <= PH.E_CUT_G and PH.E_G_H <= Eg[-1]
    lines.append("    表的范围 %.4g – %.4g eV，要盖住截断 %.4g 和最高 γ %.5g  %s"
                 % (Eg[0], Eg[-1], PH.E_CUT_G, PH.E_G_H, "✅" if c1 else "❌"))
    ok = ok and c1
    for E in (0.97 * Eg[0], 1.09 * Eg[-1]):
        try:
            PH.kn_TE_tab(pho, E)
            fired = False
        except AssertionError:
            fired = True
        lines.append("    表外 E = %.4g eV：越界 assert %s" % (E, "✅ 响了" if fired else "❌ 没响（会悄悄外推出一个错的数）"))
        ok = ok and fired
    if verbose:
        print("\n".join(lines))
    return ok


# ---------------- 自检：故意写错的版本（纪律 4） ----------------
def _te_template(E, me=ME_C2_REF, one_minus=True, num_k=True, thomson=False):
    """照 kn_T_E 的写法复刻一个，开关打开 = 故意写错。"""
    a = E / me
    mu = np.linspace(-1, 1, 2001)
    k = 1 / (1 + a * (1 - mu))
    f = (1 + mu**2) if thomson else k**2 * (k + 1 / k - 1 + mu**2)
    num = f * k if num_k else f
    trap = lambda y: (y[0] / 2 + y[1:-1].sum() + y[-1] / 2) * (2 / 2000)   # 梯形，和 photon.trap 同一种写法
    r = trap(num) / trap(f)
    return 1 - r if one_minus else r


WRONG_KT1 = [
    ("S1 少了「1 −」（算成 ⟨E′/E⟩）", lambda E: _te_template(E, one_minus=False)),
    ("S2 分子漏乘 k",                 lambda E: _te_template(E, num_k=False)),
    ("S3 mₑc² 用了 MeV（0.511）",     lambda E: _te_template(E, me=0.511)),
    ("S4 f 只用 Thomson (1 + μ²)",    lambda E: _te_template(E, thomson=True)),
]


def _fake_pho(E, TE):
    return {"KN": {"E": np.asarray(E, dtype=float), "TE": np.asarray(TE, dtype=float)}}


def selftest(pho):
    all_ok = True
    for name, fn in WRONG_KT1:
        fired = not kt1(fn, verbose=False)
        print("  %-30s %s" % (name, "✅ 响了" if fired else "❌ 没响"))
        all_ok = all_ok and fired
    quiet = kt1(lambda E: _te_template(E), verbose=False)
    print("  %-30s %s" % ("对照：KT-1 模板 + 正确公式", "✅ 不响" if quiet else "❌ 误响"))
    all_ok = all_ok and quiet

    Eg = pho["KN"]["E"]
    TE = pho["KN"]["TE"]
    shifted = _fake_pho(Eg, np.append(TE[1:], TE[-1]))                      # TE[i] 存成 E[i+1] 的值
    E10 = np.geomspace(Eg[0], Eg[-1], 10)
    coarse = _fake_pho(E10, [PH.kn_T_E(e) for e in E10])                    # 网格太粗：只有 10 点
    for name, fp in (("S5 表错位一格", shifted), ("S6 表只有 10 点", coarse)):
        fired = not kt2(fp, verbose=False)
        print("  %-30s %s" % (name, "✅ 响了" if fired else "❌ 没响"))
        all_ok = all_ok and fired
    return all_ok


if __name__ == "__main__":
    t0 = time.time()
    pho = PH.load_photon()
    print("[KT-1] 现算 kn_T_E vs 闭式解 σ_tr / σ_KN")
    ok1 = kt1()
    print("[KT-2] 表 kn_TE_tab vs 现算 kn_T_E")
    ok2 = kt2(pho)
    print("[KT-3] 表的覆盖与越界")
    ok3 = kt3(pho)
    print("自检：6 个故意写错的版本必须响")
    ok_s = selftest(pho)
    ok = ok1 and ok2 and ok3 and ok_s
    print("\n  → %s    用时 %.1f 秒" % ("✅ KT 全过" if ok else "❌ KT 没全过", time.time() - t0))
    assert ok, "KT 没全过：看上面哪一行是 ❌（KT-1 ❌ = kn_T_E 公式；KT-2 ❌ = 表；KT-3 ❌ = 表的范围或越界 assert；自检 ❌ = 判据不够灵）"
