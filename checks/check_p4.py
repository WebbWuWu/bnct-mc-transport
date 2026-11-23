# -*- coding: utf-8 -*-
# ============================================================
#  check_p4.py  --  M4 判据 P4：同一次模拟里，径迹长度估计器 Dg_t vs 碰撞估计器 Dg_c
#  判据归 Claude。2026-11-23 写。运行：Windows Ctrl+F5，要读真库，约 4–6 分钟。
#  ⚠️ 写的时候 Claude 那边没有 h5py 和真库，没在真程序上跑过（查过语法，查过「改坏」要替换的那两段原文
#     在 photon.py 里各正好出现一次）。第一次跑的完整输出发给 Claude；它自己报错（不是 assert 响）多半是 Claude 的锅。
#
#  比什么（同一次模拟里两条独立的路，逐条历史配对）：
#    路 A（碰撞估计）  Dg_c：每次碰撞按抽到的结果记（光电 E、康普顿 E − E′、电子对 E − 2mₑc²、截断剩余 E）
#    路 B（径迹长度）  Dg_t：每段飞行记 wt × 路程 × μ_tr × E（⟨T⟩/E 用和抽样同一个纯 KN 算）
#    两个估的都是「这条历史的光子在板里交给带电粒子的能量」，期望严格相等（纯 KN 对两边一样，近似 #32 不影响）。
#    每条历史 d = Σ层 Dg_c − Σ层 Dg_t，z = mean(d)/SE(d)，门槛 |z| ≤ 3（2026-11-23 写文件时定，没看过任何数据）。
#    配对比：两组来自同一批光子、正相关，配对后涨落抵掉一部分，比分别算 SE 更灵。
#  顺带核一笔账：finish 交出来的 Dg_c_sum / Dg_t_sum 30 层加起来 == 逐条历史加起来（抓 end_history、finish 接线错）。
#  逐层的比值也打印出来看，不判（R = SE/均值 > 0.2 的层不打 z，纪律 9）。
#
#  自检（纪律 4：在故意写错的版本上试响；每个用 6000 条）：
#    W3  出板那一段不记（track_layer 挪到判出板之后才调）              → 必须响
#    W4  康普顿交出比例用康普顿边 1 − 1/(1 + 2α)（最大值当平均，11/22 预测错过）→ 必须响
#    W1  光电重复计数（Dg_t 里光电项多算一遍，期望上 ≡ 吸收时又记一次）  → 只报大小、不要求响：
#        MeV 光子在水里光电份额很小，偏差可能比分辨力还小 —— 照实打印「偏多少、分辨力多少」
#    W2  截断剩余能量 Dg_t 漏记                                         → 只报大小、不要求响（每个被截断的光子 ≤ 1 keV，题库 C28）
#  抓不到什么（照实写）：
#    · 两组共用的东西一起错：截面（归 μ_tr vs NIST、OpenMC 对拍）、光子的产生（归 G5）、纯 KN 本身（近似 #32）
#    · 比分辨力小的偏差（每次打印 3 × SE ÷ 均值）
# ============================================================
import os, sys
_MY_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # my_code 根目录
for _sub in ("whiteboard", "lessons", "openmc_compare", "checks", ""):
    _p = os.path.join(_MY_CODE, _sub) if _sub else _MY_CODE
    if _p not in sys.path:
        sys.path.insert(0, _p)

import math
import time
import inspect

import transport
import photon as PH
from water import build_water
from dose import kur_table, KEEP_FACTOR

Z_GATE = 3.0
D_SLAB, NL = 15.0, 30
N_MAIN, SEED_MAIN = 20000, 40213
N_SELF = 6000
SEED_SELF = {"W1": 51829, "W2": 62417, "W3": 73901, "W4": 84253}

ORIG_TRACK = PH.track_photon
SRC_TRACK = inspect.getsource(ORIG_TRACK)
ORIG_MU_TR = PH.mu_tr
ORIG_KN_TAB = PH.kn_TE_tab


def run_capture(N, seed, nuclides, order, src):
    """跑一次 run_slab(gamma=True, 隐式)，返回每条历史的 (Σ层 Dg_c, Σ层 Dg_t) 两个列表和 run_slab 的结果。"""
    c_list, t_list = [], []
    orig_end = transport.end_history

    def end(tal, hist, nuclides_, NL_):
        c_list.append(sum(hist["Dg_c"]))
        t_list.append(sum(hist["Dg_t"]))
        return orig_end(tal, hist, nuclides_, NL_)

    transport.end_history = end
    try:
        res = transport.run_slab(N, seed, D_SLAB, nuclides, None, order, z_src=0.0, implicit=True, NL=NL,
                                 src=src, gamma=True)
    finally:
        transport.end_history = orig_end
    assert len(c_list) == N, "应该收到 %d 条历史，实际 %d 条 —— end_history 没被每条历史调一次？" % (N, len(c_list))
    return c_list, t_list, res


def paired(c_list, t_list):
    n = len(c_list)
    d = [c - t for c, t in zip(c_list, t_list)]
    md = sum(d) / n
    var = sum((x - md) ** 2 for x in d) / (n - 1)
    se = math.sqrt(var / n)
    mc = sum(c_list) / n
    mt = sum(t_list) / n
    z = md / se if se > 0 else float("inf")
    return mc, mt, md, se, z


def report(tag, c_list, t_list):
    mc, mt, md, se, z = paired(c_list, t_list)
    print("  %s：Dg_c %.1f eV/历史 · Dg_t %.1f eV/历史 · c ÷ t = %.4f" % (tag, mc, mt, mc / mt))
    print("  %s  配对差 %+.2f ± %.2f eV → z = %+.2f；偏 %+.2f%%，分辨力（3 × SE ÷ Dg_t）= %.1f%%"
          % (" " * len(tag), md, se, z, 100 * md / mt, 300 * se / mt))
    return z


def layer_table(res, N):
    """逐层看一眼（不判）。"""
    print("  逐层（不判；R > 0.2 的层不打 z）：层  Dg_c        Dg_t        c÷t     z")
    for j in (0, 2, 5, 10, 15, 20, 25, NL - 1):
        out = []
        for key in ("Dg_c", "Dg_t"):
            s, q = res[key + "_sum"][j], res[key + "_sq"][j]
            m = s / N
            se = math.sqrt(max(q / N - m * m, 0.0) / (N - 1))
            out.append((m, se))
        (mc, sc), (mt, st) = out
        if mt <= 0.0 or mc <= 0.0:
            print("        %2d  %10.3f  %10.3f   （有一组是 0，跳过）" % (j, mc, mt))
            continue
        R = max(sc / mc, st / mt)
        if R > 0.2:
            print("        %2d  %10.3f  %10.3f  %6.3f   R = %.2f > 0.2，不算 z" % (j, mc, mt, mc / mt, R))
        else:
            z = (mc - mt) / math.sqrt(sc * sc + st * st)
            print("        %2d  %10.3f  %10.3f  %6.3f  %+5.2f" % (j, mc, mt, mc / mt, z))


def spoil_track(old, new):
    """在 track_photon 的源码里把 old 换成 new，重新定义到 photon 模块里（run_photons 会调到改坏的那个）。"""
    assert SRC_TRACK.count(old) == 1, "改坏失败：track_photon 里「%s」出现了 %d 次（应该正好 1 次）—— 你改过这一行？告诉 Claude" % (old, SRC_TRACK.count(old))
    exec(compile(SRC_TRACK.replace(old, new), "<改坏的 track_photon>", "exec"), PH.__dict__)
    assert PH.track_photon is not ORIG_TRACK, "改坏没生效"


def restore():
    PH.track_photon = ORIG_TRACK
    PH.mu_tr = ORIG_MU_TR
    PH.kn_TE_tab = ORIG_KN_TAB


def self_test(name, desc, must_ring, setup, nuclides, order, src):
    print("自检 %s：%s，%d 条，%s" % (name, desc, N_SELF, "必须响" if must_ring else "只报大小"))
    setup()
    try:
        c, t, _ = run_capture(N_SELF, SEED_SELF[name], nuclides, order, src)
    finally:
        restore()
    z = report(name, c, t)
    if must_ring:
        rang = abs(z) > Z_GATE
        print("  %s\n" % ("✅ 响了" if rang else "❌ 没响 —— P4 抓不住这个错，告诉 Claude"))
        return rang
    print("  （不判：%s）\n" % ("这么大的偏差 6000 条就分辨出来了" if abs(z) > Z_GATE else "比分辨力小，P4 抓不住这一类，照实记进判据说明"))
    return True


if __name__ == "__main__":
    t0 = time.time()
    edges, shares = kur_table()
    nuclides, order = build_water(KEEP_FACTOR * edges[-1])
    src = (edges, shares)
    ok = True

    print("P4 主判据：水模体 + KUR 源，隐式俘获，%d 条" % N_MAIN)
    c, t, res = run_capture(N_MAIN, SEED_MAIN, nuclides, order, src)
    if sum(t) == 0.0 or sum(c) == 0.0:
        print("  ❌ 有一组一点都没记（Dg_c 总和 %.3g，Dg_t 总和 %.3g）—— 光子没跑起来？下面不判了" % (sum(c), sum(t)))
        sys.exit(1)
    z = report("主", c, t)
    good = abs(z) <= Z_GATE
    ok = ok and good
    print("  %s |z| ≤ %.0f" % ("✅" if good else "❌", Z_GATE))

    acc_ok = True
    for key, lst in (("Dg_c", c), ("Dg_t", t)):
        tot_tal = sum(res[key + "_sum"])
        tot_his = sum(lst)
        same = abs(tot_tal - tot_his) <= 1e-9 * abs(tot_his)
        acc_ok = acc_ok and same
        print("  对账 %s：finish 交出的 30 层之和 %.6e · 逐条历史之和 %.6e  %s" % (key, tot_tal, tot_his, "✅" if same else "❌"))
    ok = ok and acc_ok
    layer_table(res, N_MAIN)
    print()

    tests = [
        ("W3", "出板那一段不记（track_layer 挪到判出板之后）", True,
         lambda: spoil_track('track_layer(hist["Dg_t"]', '(0 <= z1 <= d) and track_layer(hist["Dg_t"]')),
        ("W4", "康普顿交出比例用康普顿边（最大值当平均）", True,
         lambda: setattr(PH, "kn_TE_tab", lambda pho, E: 1.0 - 1.0 / (1.0 + 2.0 * E / PH.ME_C2))),
        ("W1", "光电重复计数（Dg_t 的光电项多算一遍）", False,
         lambda: setattr(PH, "mu_tr", lambda sig, E, pho: ORIG_MU_TR(sig, E, pho) + sig["pe"])),
        ("W2", "截断剩余能量 Dg_t 漏记", False,
         lambda: spoil_track('hist["Dg_t"][k]+=g["wt"]*g["E"]', 'pass')),
    ]
    for name, desc, must, setup in tests:
        ok = self_test(name, desc, must, setup, nuclides, order, src) and ok

    print("用时 %.1f 分钟" % ((time.time() - t0) / 60))
    assert ok, "P4 没全过：看上面哪一行是 ❌（主判据 ❌ = 两组记的不是同一个量；对账 ❌ = end_history / finish 接线；自检 ❌ = 判据不够灵）"
    print("P4 全部通过")
