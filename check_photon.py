# -*- coding: utf-8 -*-
# ============================================================
#  check_photon.py  --  M4 判据 P1：光子截面读得对不对、水合成得对不对
#  Claude 2026-11-08 写（判据归 Claude，成品直接跑，不用默写）
#
#  用法：在 VS Code 里打开本文件，Ctrl+F5（Windows 的 openmc (3.10.20) 环境）。
#        它自己搭水模体（lesson11_water.build_water，要读中子截面，十几秒）、读光子库（photon.load_photon），
#        先跑 check_P1（正式判），再跑 check_P1_selftest（5 种故意写错的必须响、响对的那一条，纪律 4）。
#
#  两条路（纪律 2）：
#    你的路   = photon.photon_sigma(pho, nuclides, E)：ENDF/B-VIII.0 光原子库（EPDL97 系），按元素读、按个数密度合成
#    判据的路 = data\photon\water_xcom.txt：NIST XCOM 直接给的「水」这个化合物的 μ/ρ（另一套评价、另一种合成）
#    先确认比的是同一个量：XCOM 第 6 列 = 总衰减（含相干）cm²/g；你的 Σ_tot（/cm）÷ ρ = 0.998 g/cm³ 就是它
#  3 条：
#    [P1-1] 6 个能量点（0.03、0.1、0.4778、1、2.223、3 MeV，都是 XCOM 表里原有的点，不用插值）：
#           Σ_tot / ρ vs XCOM 第 6 列，|相对差| < 2%                     抓：漏乘 1e-24、= 写成 +=、个数密度 / 元素对错
#    [P1-2] 0.4778 MeV 查电子对：不报错、正好 = 0；2.223 MeV 电子对占比 vs XCOM（1.176%），|差| < 0.15 个百分点
#                                                                        抓：log(0) 没处理（报错）、阈值附近插错
#    [P1-3] 0.03、0.1、1 MeV：四项各自占 Σ_tot 的份额 vs XCOM，|差| < 2 个百分点
#                                                                        抓：反应名字对错位（总量不变，P1-1 抓不到）
#  门槛在看数据之前定（11/8 18:20）：两套评价本来差约 1%，门槛取 2%；份额用百分点（小份额时相对差没意义）
#  盲区：只查 6 个能量点；1 keV 以下、10 MeV 以上不查（光子用不到）；B 的份额太小（15 ppm），B 读错这里看不出
# ============================================================

import os
import math
import numpy as np

from photon import load_photon, photon_sigma, loglog0, NUC2ELEM
from lesson11_water import build_water
from lesson08_dose import kur_table, KEEP_FACTOR

HERE = os.path.dirname(os.path.abspath(__file__))
XCOM_PATH = os.path.join(HERE, "data", "photon", "water_xcom.txt")
RHO = 0.998                                       # g/cm³，和 lesson11_water.rho_w 一致
E_P1 = [0.03, 0.1, 0.4778, 1.0, 2.223, 3.0]       # MeV，XCOM 表里原有的点
E_P3 = [0.03, 0.1, 1.0]                           # MeV
TOL_TOT = 0.02                                    # P1-1：相对差
TOL_PAIR = 0.0015                                 # P1-2：0.15 个百分点
TOL_SHARE = 0.02                                  # P1-3：2 个百分点
KEYS = ["coh", "incoh", "pe", "pair"]


def read_xcom():
    """读 XCOM 表 → {能量 MeV: {"coh","incoh","pe","pair","tot"}}，单位 cm²/g。"""
    d = np.loadtxt(XCOM_PATH, comments="#", encoding="utf-8")
    out = {}
    for row in d:
        out[round(float(row[0]), 6)] = {"coh": row[1], "incoh": row[2], "pe": row[3],
                                        "pair": row[4] + row[5], "tot": row[6]}
    return out


def _xc(xcom, E_MeV):
    key = round(E_MeV, 6)
    assert key in xcom, "XCOM 表里没有 %g MeV 这个点（判据自己的错）" % E_MeV
    return xcom[key]


def check_P1(sigma_func, pho, nuclides, xcom, verbose=True):
    """返回不过的条目列表（空 = 全过）。sigma_func 签名同 photon_sigma。"""
    bad = []
    # ---------- P1-1 ----------
    worst = 0.0
    for E in E_P1:
        try:
            s = sigma_func(pho, nuclides, E * 1e6)
        except Exception as e:
            bad.append("[P1-1] %g MeV 报错：%r" % (E, e))
            continue
        mine = s["tot"] / RHO
        ref = _xc(xcom, E)["tot"]
        r = mine / ref - 1
        worst = max(worst, abs(r))
        if verbose:
            print("  [P1-1] %6.4g MeV · 你 %.5g · XCOM %.5g cm²/g · 差 %+.2f%%" % (E, mine, ref, 100 * r))
        if abs(r) >= TOL_TOT:
            bad.append("[P1-1] %g MeV：μ/ρ 你 %.5g vs XCOM %.5g，差 %+.1f%%（门槛 ±2%%）" % (E, mine, ref, 100 * r))
    # ---------- P1-2 ----------
    try:
        s = sigma_func(pho, nuclides, 0.4778e6)
        if s["pair"] != 0.0:
            bad.append("[P1-2] 0.4778 MeV 电子对 = %.3g /cm，应正好是 0（阈值 1.022 MeV 以下）" % s["pair"])
        s2 = sigma_func(pho, nuclides, 2.223e6)
        x2 = _xc(xcom, 2.223)
        mine_sh = s2["pair"] / s2["tot"]
        ref_sh = x2["pair"] / x2["tot"]
        if verbose:
            print("  [P1-2] 0.4778 MeV 电子对 = %g · 2.223 MeV 电子对占比 你 %.3f%% · XCOM %.3f%%"
                  % (s["pair"], 100 * mine_sh, 100 * ref_sh))
        if abs(mine_sh - ref_sh) >= TOL_PAIR:
            bad.append("[P1-2] 2.223 MeV 电子对占比 你 %.3f%% vs XCOM %.3f%%（门槛 ±0.15 个百分点）"
                       % (100 * mine_sh, 100 * ref_sh))
    except Exception as e:
        bad.append("[P1-2] 查电子对报错：%r（截面 = 0 时 log(0)？）" % (e,))
    # ---------- P1-3 ----------
    for E in E_P3:
        try:
            s = sigma_func(pho, nuclides, E * 1e6)
        except Exception as e:
            bad.append("[P1-3] %g MeV 报错：%r" % (E, e))
            continue
        x = _xc(xcom, E)
        line = []
        for k in KEYS:
            mine = s[k] / s["tot"]
            ref = x[k] / x["tot"]
            line.append("%s %.1f/%.1f" % (k, 100 * mine, 100 * ref))
            if abs(mine - ref) >= TOL_SHARE:
                bad.append("[P1-3] %g MeV：%s 占比 你 %.1f%% vs XCOM %.1f%%（门槛 ±2 个百分点）"
                           % (E, k, 100 * mine, 100 * ref))
        if verbose:
            print("  [P1-3] %5.3g MeV 份额（你/XCOM，%%）：%s" % (E, " · ".join(line)))
    return bad


# ---------------- 故意写错的 5 个版本（纪律 4） ----------------
def _loglog_no0(x, y, e):
    """错：不判 0，直接取 log。"""
    lo = np.searchsorted(x, e, side="right") - 1
    frac = (math.log(e) - math.log(x[lo])) / (math.log(x[lo + 1]) - math.log(x[lo]))
    return math.exp(math.log(y[lo]) + frac * (math.log(y[lo + 1]) - math.log(y[lo])))


def _make_wrong(kind):
    def f(photon, nuclides, E):
        sig = {k: 0.0 for k in KEYS}
        for nuc in nuclides:
            elem = NUC2ELEM[nuc["name"]]
            if kind == "swapHO":
                elem = {"H": "O", "O": "H", "B": "B"}[elem]
            for k in KEYS:
                kk = k
                if kind == "swap_coh_pe":
                    kk = {"coh": "pe", "pe": "coh"}.get(k, k)
                interp = _loglog_no0 if kind == "no0" else loglog0
                v = interp(photon[elem]["E"], photon[elem][kk], E) * nuc["dens"]
                v = v if kind == "no_barn" else v * 1e-24
                if kind == "assign":
                    sig[k] = v
                else:
                    sig[k] += v
        sig["tot"] = sum(sig[k] for k in KEYS)
        return sig
    return f


WRONG = [("漏乘 1e-24", "no_barn", "[P1-1]"),
         ("+= 写成 =（只剩最后一个核素 B10）", "assign", "[P1-1]"),
         ("H、O 元素对调", "swapHO", "[P1-1]"),
         ("相干和光电名字对调", "swap_coh_pe", "[P1-3]"),
         ("截面为 0 时不判、直接 log", "no0", "[P1-2]")]


def check_P1_selftest(pho, nuclides, xcom):
    ok = True
    for name, kind, expect in WRONG:
        bad = check_P1(_make_wrong(kind), pho, nuclides, xcom, verbose=False)
        hit = [b for b in bad if b.startswith(expect)]
        if hit:
            print("自检 ✅ %-22s → %s 响了" % (name, expect))
        else:
            ok = False
            print("自检 ❌ %-22s → %s 没响（实际：%s）" % (name, expect, bad[:2]))
    assert ok, "判据 P1 自检没全过：有的错没抓到 —— 这是判据的问题，发给 Claude"
    print("判据 P1 自检全过：5 种错都响了，而且响的是对的那一条")


if __name__ == "__main__":
    edges, shares = kur_table()
    nuclides, order = build_water(KEEP_FACTOR * edges[-1])
    pho = load_photon()
    xcom = read_xcom()
    print("水模体 + 光子库读好；核素：%s" % [n["name"] for n in nuclides])
    bad = check_P1(photon_sigma, pho, nuclides, xcom)
    for b in bad:
        print("❌", b)
    assert not bad, "判据 P1 没过（%d 处，见上）" % len(bad)
    print("判据 P1 过")
    check_P1_selftest(pho, nuclides, xcom)
