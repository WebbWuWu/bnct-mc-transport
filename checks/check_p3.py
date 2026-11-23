# -*- coding: utf-8 -*-
# ============================================================
#  check_p3.py  --  M4 判据 P3：单能 2.2233 MeV 光子垂直入射 15 cm 水，一次不撞就穿过的份额 vs 闭式解
#  判据归 Claude。2026-11-23 写。运行：Windows Ctrl+F5，要读光子库，约半分钟到一分钟。
#  ⚠️ 写的时候 Claude 那边没有 h5py 和真库，没在真程序上跑过（只查过语法）。
#     第一次跑的完整输出发给 Claude；如果它自己报错（不是 assert 响），多半是 Claude 的锅。
#
#  比什么（两条独立的路）：
#    路 A（模拟）  N 个 2.2233 MeV 光子从 z = 0 沿 +z 射进去，交给你的 track_photon 跑。
#                  怎么知道「一次都没撞」：数这个光子用了几个随机数（rng.n）。
#                  只用了 1 个 = 只抽了一次飞行距离就出板了；撞过的至少还要抽「哪种反应」那一个。
#    路 B（闭式解）p = exp(−Σ_tot · 15 cm)，Σ_tot 取 photon_sigma 在 2.2233 MeV 的值。
#    数到的份额是二项分布，SE = √(p(1 − p)/N)；z = (份额 − p)/SE，门槛 |z| ≤ 3（写文件时定，没看过数据）。
#  再加一条量级检查（第三条路）：Σ_tot 和 NIST XCOM 水的 μ/ρ × 0.998 比，相对差 < 2%（同样先定）。
#
#  抓什么：
#    · 自由程抽样写错（−ln ξ 乘 Σ 而不是除、用了 Σ_pe 而不是 Σ_tot、ξ 用了两次）
#    · 出板判断写错（z1 > d 漏判、出板后还继续撞）
#    · 截面量级错（漏乘 1e-24、密度错）—— 这条靠 NIST 那一行
#  抓不到什么（照实写）：
#    · 撞了以后的事（反应怎么选、角度、能量）—— 归 P2、P4 和 OpenMC 对拍
#    · Σ_tot 里四项怎么分 —— 这里只看总和
#  自检（纪律 4）：把 Σ_tot 改大 10%（「截面偏大」那类错）、闭式解还用原来的 Σ → 必须响
# ============================================================
import os, sys
_MY_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # my_code 根目录
for _sub in ("whiteboard", "lessons", "openmc_compare", "checks", ""):
    _p = os.path.join(_MY_CODE, _sub) if _sub else _MY_CODE
    if _p not in sys.path:
        sys.path.insert(0, _p)

import math
import time

import photon as PH
from rng import lcg
from water import build_water
from dose import kur_table, KEEP_FACTOR

D_SLAB, NL = 15.0, 30
E0 = PH.E_G_H                      # 2.2233 MeV
N_MAIN, SEED_MAIN = 20000, 91237
N_SELF, SEED_SELF = 8000, 23371
Z_GATE = 3.0
NIST_GATE = 0.02
RHO_W = 0.998                      # water.py 里的 rho_w

# NIST XCOM 水 μ/ρ（总，含相干），cm²/g；2.2233 MeV 在 2 和 3 MeV 之间按 log-log 插
NIST_E = (2.0e6, 3.0e6)
NIST_MU = (4.942e-2, 3.969e-2)


def nist_mu(E):
    f = math.log(E / NIST_E[0]) / math.log(NIST_E[1] / NIST_E[0])
    return math.exp(math.log(NIST_MU[0]) + f * (math.log(NIST_MU[1]) - math.log(NIST_MU[0])))


def count_through(N, seed, pho, nuclides, E, scale_tot=1.0):
    """N 个光子垂直入射，返回一次没撞就穿过的个数。scale_tot ≠ 1：自检，把 Σ_tot 放大。"""
    rng = lcg(seed)
    dz = D_SLAB / NL
    orig_sigma = PH.photon_sigma

    def sigma_spoiled(photon, nuclides_, E_):
        sig = orig_sigma(photon, nuclides_, E_)
        sig["tot"] *= scale_tot
        return sig

    if scale_tot != 1.0:
        PH.photon_sigma = sigma_spoiled
    n_through = 0
    try:
        for i in range(N):
            g = {"z": 0.0, "u": 0.0, "v": 0.0, "w": 1.0, "E": E, "wt": 1.0}
            hist = {"Dg_c": [0.0] * NL, "Dg_t": [0.0] * NL, "pstack": [], "gsites": []}
            n0 = rng.n
            PH.track_photon(g, hist, rng, pho, nuclides, D_SLAB, dz, NL)
            used = rng.n - n0
            assert used >= 1, "第 %d 个光子一个随机数都没用 —— 连飞行距离都没抽？" % i
            if used == 1:
                n_through += 1
    finally:
        PH.photon_sigma = orig_sigma
    return n_through


def judge(tag, n, N, p):
    frac = n / N
    se = math.sqrt(p * (1.0 - p) / N)
    z = (frac - p) / se
    print("  %s：穿过 %d / %d = %.5f · 闭式解 %.5f · 差 %+.5f ± %.5f → z = %+.2f" % (tag, n, N, frac, p, frac - p, se, z))
    return z


if __name__ == "__main__":
    t0 = time.time()
    edges, shares = kur_table()
    nuclides, order = build_water(KEEP_FACTOR * edges[-1])
    pho = PH.load_photon()
    ok = True

    sig_tot = PH.photon_sigma(pho, nuclides, E0)["tot"]
    p = math.exp(-sig_tot * D_SLAB)
    ref = nist_mu(E0) * RHO_W
    rel = sig_tot / ref - 1.0
    good = abs(rel) < NIST_GATE
    ok = ok and good
    print("量级：Σ_tot(%.4f MeV) 程序 %.5f /cm · NIST × 0.998 = %.5f /cm · 相对差 %+.2f%%  %s（门槛 %.0f%%）"
          % (E0 / 1e6, sig_tot, ref, 100 * rel, "✅" if good else "❌", 100 * NIST_GATE))
    print("闭式解：p = exp(−%.5f × %.0f) = %.5f\n" % (sig_tot, D_SLAB, p))

    print("P3 主判据：%d 个光子" % N_MAIN)
    n = count_through(N_MAIN, SEED_MAIN, pho, nuclides, E0)
    z = judge("主", n, N_MAIN, p)
    good = abs(z) <= Z_GATE
    ok = ok and good
    print("  %s |z| ≤ %.0f\n" % ("✅" if good else "❌", Z_GATE))

    print("自检：Σ_tot × 1.10（闭式解还用原来的 Σ），%d 个，必须响" % N_SELF)
    n2 = count_through(N_SELF, SEED_SELF, pho, nuclides, E0, scale_tot=1.10)
    z2 = judge("改坏", n2, N_SELF, p)
    rang = abs(z2) > Z_GATE
    ok = ok and rang
    print("  %s\n" % ("✅ 响了" if rang else "❌ 没响 —— 条数太少分辨不出 10%，告诉 Claude"))

    print("用时 %.1f 分钟" % ((time.time() - t0) / 60))
    assert ok, "P3 没全过：看上面哪一行是 ❌（量级 ❌ = 截面单位 / 密度；主判据 ❌ = 自由程抽样或出板判断；自检 ❌ = 判据不够灵）"
    print("P3 全部通过")
