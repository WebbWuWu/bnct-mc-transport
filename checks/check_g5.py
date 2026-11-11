# -*- coding: utf-8 -*-
# ============================================================
#  check_g5.py  --  M4 判据 G5：生出来的光子总权重，对不对得上中子那边「该出多少 γ」
#  判据归 Claude。2026-11-11 写（11/9 用户定「G5 加」，一直欠着）。运行：Windows Ctrl+F5，要读真库，约 3–4 分钟。
#  ⚠️ 写的时候 Claude 那边没有 h5py 和真库，本文件没在真程序上跑过（只查过语法、用假模块试过流程）。
#     第一次跑的完整输出发给 Claude；如果它自己报错（不是 assert 响），多半是 Claude 的锅。
#
#  比什么（同一次模拟里的两条独立的路，逐条历史配对）：
#    路 A（碰撞估计）  你的 record_site + emit_photons：每次碰撞按 y = (Σ_γH + Σ_801B)/Σ_t 生 ⌊y + ξ⌋ 个光子，
#                      每个带中子的权重 → 一条历史生出来的光子权重之和 a
#    路 B（径迹长度）  M2 起就有的 γ 计数：fly 里每段径迹记 wt × ℓ × (Σ_γH + Σ_801B) → 一条历史 30 层加起来 b
#    两个都在估「这条历史平均放出多少个俘获 γ（按权重算）」，期望相同。
#    每条历史算 d = a − b，N 条的均值应该是 0：z = mean(d) / SE(d)，门槛 |z| ≤ 3（🔴 2026-11-11 写文件时定，没看过任何数据）。
#    用配对差而不是分别算 SE：a 和 b 来自同一条历史、正相关，配对以后涨落互相抵掉一部分，分辨力更高。
#
#  抓什么：
#    · 光子权重没继承中子权重（题库 C2：每个光子权重设成 1）
#    · record_site 记的是打折 / 轮盘赌之后的权重（应该是碰撞前的）
#    · y 的分子分母写错、漏了硼那条线、⌊y + ξ⌋ 写成四舍五入
#    · run_photons 漏生 / 重复生（比如产生点列表没清空，下一条历史又生一遍）
#  抓不到什么（照实写）：
#    · 2.2233 / 0.4778 MeV 两条线的份额、光子方向、位置 z —— 这些归 G1、G4
#    · 路 A、路 B 共用的 Σ_γH、Σ_801B 本身错了（两边一起错）—— 这个归 M2 的对拍
#    · 小于分辨力的偏差：分辨力每次运行都打印（3 × SE ÷ 均值），2 万条时大约 2–3%
#
#  自检（判据要在故意写错的版本上试响）：
#    S1  把每个新生光子的权重改成 1.0（题库 C2 那个错），隐式俘获下重跑 3000 条 → 必须响
#    S2  把路 A 整体乘 0.95（等于「少生 5% 的光子」），用主运行的数直接算 → 必须响；不响说明条数太少、G5 分辨不出 5%
#  模拟模式（implicit=False）也跑一次：那里中子权重恒为 1，S1 那个错在模拟模式下看不出来 —— 所以主判据放在隐式模式。
# ============================================================
# ---- 路径：本文件在 my_code 的子文件夹里，把 my_code 根目录和各子文件夹加进搜索路径（和 check_emit.py 同一段）----
import os, sys
_MY_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # my_code 根目录
for _sub in ("whiteboard", "lessons", "openmc_compare", "checks", ""):
    _p = os.path.join(_MY_CODE, _sub) if _sub else _MY_CODE
    if _p not in sys.path:
        sys.path.insert(0, _p)

import math
import time

import transport
import photon as PH
from water import build_water
from dose import kur_table, KEEP_FACTOR

Z_GATE = 3.0
D_SLAB, NL = 15.0, 30
N_MAIN, SEED_MAIN = 20000, 60719        # 主判据：隐式
N_ANALOG, SEED_ANALOG = 8000, 35281     # 模拟模式
N_SELF, SEED_SELF = 3000, 77317         # 自检 S1


def run_pair(N, seed, implicit, nuclides, order, src, spoil_weight=False):
    """跑一次 run_slab(gamma=True)，返回每条历史的 (a, b) 两个列表。
    a = 这条历史 emit_photons 压进栈的光子权重之和；b = 这条历史 30 层 gH_h + gB_h 之和。
    spoil_weight=True：自检用，把每个新生光子的权重改成 1.0。"""
    a_list, b_list = [], []
    cur = [0.0]
    orig_emit = PH.emit_photons
    orig_end = transport.end_history

    def emit(site, hist, rng):
        before = len(hist["pstack"])
        n = orig_emit(site, hist, rng)
        new = hist["pstack"][before:]
        assert len(new) == n, "emit_photons 返回 %d，但栈里多了 %d 个" % (n, len(new))
        if spoil_weight:
            for g in new:
                g["wt"] = 1.0
        cur[0] += sum(g["wt"] for g in new)
        return n

    def end(tal, hist, nuclides_, NL_):
        a_list.append(cur[0])
        b_list.append(sum(hist["gH_h"]) + sum(hist["gB_h"]))
        cur[0] = 0.0
        return orig_end(tal, hist, nuclides_, NL_)

    PH.emit_photons = emit
    transport.end_history = end
    try:
        transport.run_slab(N, seed, D_SLAB, nuclides, None, order, z_src=0.0, implicit=implicit, NL=NL,
                           src=src, gamma=True)
    finally:
        PH.emit_photons = orig_emit
        transport.end_history = orig_end
    assert len(a_list) == N, "应该收到 %d 条历史，实际 %d 条 —— end_history 没被每条历史调一次？" % (N, len(a_list))
    return a_list, b_list


def paired(a_list, b_list, scale_a=1.0):
    """配对差 d = scale_a·a − b 的均值、标准误、z；顺带给 a、b 各自的均值。"""
    n = len(a_list)
    d = [scale_a * a - b for a, b in zip(a_list, b_list)]
    md = sum(d) / n
    var = sum((x - md) ** 2 for x in d) / (n - 1)
    se = math.sqrt(var / n)
    ma = scale_a * sum(a_list) / n
    mb = sum(b_list) / n
    z = md / se if se > 0 else float("inf")
    return ma, mb, md, se, z


def report(tag, a_list, b_list, scale_a=1.0):
    ma, mb, md, se, z = paired(a_list, b_list, scale_a)
    print("  %s：路 A（生出来的光子权重）%.5f / 历史 · 路 B（径迹长度 γ 计数）%.5f / 历史 · A ÷ B = %.4f" % (tag, ma, mb, ma / mb))
    print("  %s  配对差 %.5f ± %.5f → z = %+.2f；分辨力（3 × SE ÷ B）= %.1f%%" % (" " * len(tag), md, se, z, 3 * se / mb * 100))
    return z, mb


if __name__ == "__main__":
    edges, shares = kur_table()
    nuclides, order = build_water(KEEP_FACTOR * edges[-1])
    src = (edges, shares)
    ok = True
    t0 = time.time()

    print("G5 主判据：水模体 + KUR 源，隐式俘获，%d 条（约 2 分钟）" % N_MAIN)
    a, b = run_pair(N_MAIN, SEED_MAIN, True, nuclides, order, src)
    if sum(a) == 0.0:
        print("  ❌ 一个光子都没生 —— run_photons 没调 emit_photons？gamma 没传下去？下面不判了")
        sys.exit(1)
    z, mb = report("隐式", a, b)
    good = abs(z) <= Z_GATE
    ok = ok and good
    print("  %s |z| ≤ %.0f\n" % ("✅" if good else "❌", Z_GATE))

    print("自检 S2（用上面的数直接算）：路 A 整体 × 0.95，必须响")
    _, _, _, _, z2 = paired(a, b, 0.95)
    rang = abs(z2) > Z_GATE
    ok = ok and rang
    print("  %s z = %+.2f%s\n" % ("✅ 响了" if rang else "❌ 没响", z2, "" if rang else " —— %d 条分辨不出 5%%，加条数" % N_MAIN))

    print("G5 模拟模式：%d 条（中子权重恒为 1）" % N_ANALOG)
    a2, b2 = run_pair(N_ANALOG, SEED_ANALOG, False, nuclides, order, src)
    z, _ = report("模拟", a2, b2)
    good = abs(z) <= Z_GATE
    ok = ok and good
    print("  %s |z| ≤ %.0f\n" % ("✅" if good else "❌", Z_GATE))

    print("自检 S1：把每个新生光子的权重改成 1.0（题库 C2 的错），隐式，%d 条，必须响" % N_SELF)
    a3, b3 = run_pair(N_SELF, SEED_SELF, True, nuclides, order, src, spoil_weight=True)
    z3, _ = report("改坏", a3, b3)
    rang = abs(z3) > Z_GATE
    ok = ok and rang
    print("  %s\n" % ("✅ 响了" if rang else "❌ 没响 —— G5 抓不住「光子权重 = 1」，判据要重写，告诉 Claude"))

    print("用时 %.1f 分钟" % ((time.time() - t0) / 60))
    assert ok, "G5 没全过：看上面哪一行是 ❌（主判据 / 模拟模式 ❌ = 生光子的权重或个数有问题；自检 ❌ = 判据自己不够灵）"
    print("G5 全部通过")
