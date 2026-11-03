# -*- coding: utf-8 -*-
# ============================================================
#  check_elastic.py  --  M3 判据 A：load_elastic_angle 读库、切表对不对
#  Claude 2026-11-03 写（判据归 Claude，成品直接 import，不用默写）
#
#  用法：在 elastic_angle.py 的 if __name__ == "__main__": 块里
#      from check_elastic import check_A, check_A_selftest
#      angle = load_elastic_angle(LIB, NUCS)
#      check_A(angle)              # 正式判
#      check_A_selftest(angle)     # 头一次跑做一遍：故意改错 5 种，确认每种都响、响对的那一条（纪律 4）
#  运行：Windows 的 openmc (3.10.20) 环境，Ctrl+F5（只用 h5py + numpy）
#
#  两条路：
#    你的路   = 按 mu 数据集上的 offsets 切表
#    判据的路 = 不用 offsets：在 μ 那一行里找「μ = −1」的位置当每张表的起点
#               （11/1 判据 A2 实测：每张表 μ 从 −1 单调升到 +1，所以 −1 只出现在表头）
#  4 条（跑的顺序 A-1 → A-4 → A-3 → A-2）：
#    [A-1] 张数、入射能量：张数 = 能量个数 = 判据找到的 −1 个数；能量逐个相等
#    [A-4] 质心系标志：H1、O16 都是 True（下次抽样要据此从质心系换到实验室系）
#    [A-3] 物理锚：μ̄_cm（对表精确积分）和 11/1 读出来的数对上
#          O16 0.5 MeV +0.417、O16 2 MeV +0.121、H1 2 MeV ≈ 0（|μ̄| < 0.01）
#    [A-2] 逐位：每张表的 μ、pdf、cdf 三行和判据的路切出来的完全相等（np.array_equal）
#  门槛在看数据之前定：A-3 容差 0.002（11/1 打印到第 3 位）；其余都是「完全相等」
#  盲区：两条路读的是同一个库文件 —— 库本身错了抓不到（这由整体对拍 vs OpenMC 管）
# ============================================================

import os
import copy
import h5py
import numpy as np

LIB = r"E:\BNCT_GPU\lib80x_hdf5"
A_CM_WANT = {"H1": True, "O16": True}
A_MUBAR = [("O16", 5.0e5, 0.417, 0.002),        # (核素, 入射能量 eV, 11/1 读出的 μ̄_cm, 容差)
           ("O16", 2.0e6, 0.121, 0.002),
           ("H1", 2.0e6, 0.0, 0.01)]


def a_read(nuc, lib=LIB):
    """判据自己的路：读入射能量、质心系标志、mu 那一长条（三行），不读 offsets。"""
    with h5py.File(os.path.join(lib, nuc + ".h5"), "r") as f:
        rx = f[nuc]["reactions"]["reaction_002"]
        cm = bool(rx.attrs["center_of_mass"])
        ang = rx["product_0"]["distribution_0"]["angle"]
        E = ang["energy"][:]
        mu3 = ang["mu"][0:3, :]
    return np.asarray(E, float), cm, np.asarray(mu3, float)


def a_split(mu3):
    """按「μ = −1」找表头切表。"""
    starts = np.flatnonzero(mu3[0] == -1.0)
    ends = np.append(starts[1:], mu3.shape[1])
    return [mu3[:, s:e] for s, e in zip(starts, ends)]


def a_mubar(mu, pdf):
    """一张表的 μ̄：pdf 每段是直线 → μ·p 是二次 → 辛普森精确（11/1 同一个式子），再除以 ∫p。"""
    a, b = mu[:-1], mu[1:]
    pa, pb = pdf[:-1], pdf[1:]
    norm = np.sum((b - a) * (pa + pb) / 2.0)
    return float(np.sum((b - a) / 6.0 * ((2 * a + b) * pa + (a + 2 * b) * pb)) / norm)


def check_A(angle, quiet=False, reader=a_read):
    raw = {}
    for nuc in A_CM_WANT:
        assert nuc in angle, "[A-1] 字典里没有 %s（键是不是加错了引号 / 用了别的名字？）现有的键：%r" % (nuc, list(angle.keys()))
        E, cm, mu3 = reader(nuc)
        parts = a_split(mu3)
        raw[nuc] = (E, cm, parts)
        mine = angle[nuc]
        # [A-1] 抓：offsets 用错导致张数不对、切表循环少一圈 / 多一圈、"E" 存成了别的数组
        assert len(mine["tables"]) == len(E) == len(parts), \
            "[A-1] %s：你的表 %d 张，入射能量 %d 个，判据按 −1 找到 %d 张" % (nuc, len(mine["tables"]), len(E), len(parts))
        assert np.array_equal(np.asarray(mine["E"], float), E), "[A-1] %s：\"E\" 和库里的入射能量不一样" % nuc
    if not quiet:
        print("A-1 ✅ 张数、入射能量：H1 %d 张、O16 %d 张" % (len(raw["H1"][0]), len(raw["O16"][0])))

    for nuc, want in A_CM_WANT.items():
        # [A-4] 抓：cm 读错属性 / 没转 bool / 存反
        assert angle[nuc]["cm"] is want, "[A-4] %s：\"cm\" = %r，库里 center_of_mass = %r" % (nuc, angle[nuc]["cm"], want)
    if not quiet:
        print("A-4 ✅ 质心系标志：H1、O16 都是 True")

    for nuc, E_pick, want, tol in A_MUBAR:
        E = raw[nuc][0]
        j = int(np.argmin(np.abs(E - E_pick)))
        t = angle[nuc]["tables"][j]
        mb = a_mubar(np.asarray(t["mu"], float), np.asarray(t["pdf"], float))
        # [A-3] 抓：pdf / cdf 行拿反、表号错位（拿到了别的能量的表）、μ 和 pdf 不是同一张表的
        assert abs(mb - want) < tol, \
            "[A-3] %s 第 %d 张（E = %.4g eV）μ̄_cm = %+.4f，11/1 读出来是 %+.3f（容差 %.3f）" % (nuc, j, E[j], mb, want, tol)
        if not quiet:
            print("A-3 ✅ %-3s 第 %4d 张 E = %.4g eV · μ̄_cm = %+.4f（11/1：%+.3f）" % (nuc, j, E[j], mb, want))

    for nuc in A_CM_WANT:
        parts = raw[nuc][2]
        for j, (t, p) in enumerate(zip(angle[nuc]["tables"], parts)):
            for r, key in enumerate(("mu", "pdf", "cdf")):
                # [A-2] 抓：切片起止错一格、三行的顺序拆错、某一张表放错位置
                assert np.array_equal(np.asarray(t[key], float), p[r]), \
                    "[A-2] %s 第 %d 张的 \"%s\" 和判据按 −1 切出来的不一样（长度 %d vs %d）" % (nuc, j, key, len(t[key]), len(p[r]))
    if not quiet:
        print("A-2 ✅ 逐位：两个核素全部表的 μ、pdf、cdf 三行和判据的路完全相等")
        print("判据 A 过")


def check_A_selftest(angle, reader=a_read):
    """纪律 4：在故意改错的副本上跑，必须响，而且响对的那一条。不改你的 angle（全部在深拷贝上改）。"""
    def drop_last(a):                      # 切表循环少一圈
        a["O16"]["tables"].pop()

    def cm_false(a):                       # 质心系标志存错
        a["H1"]["cm"] = False

    def iso_05(a):                         # O16 0.5 MeV 那张换成各向同性（pdf 全 0.5）
        E = np.asarray(a["O16"]["E"], float)
        j = int(np.argmin(np.abs(E - 5.0e5)))
        t = a["O16"]["tables"][j]
        t["pdf"] = np.full(len(t["mu"]), 0.5)

    def swap_rows(a):                      # 一张离锚点很远的表 pdf、cdf 拆反了（A-3 看不到 → 只有 A-2 能抓）
        t = a["H1"]["tables"][7]
        t["pdf"], t["cdf"] = t["cdf"], t["pdf"]

    def shift_one(a):                      # 一张表切片终点多 +1（多吞了下一张表的第一个点）
        nxt = a["O16"]["tables"][101]
        t = a["O16"]["tables"][100]
        t["mu"] = np.append(t["mu"], nxt["mu"][0])
        t["pdf"] = np.append(t["pdf"], nxt["pdf"][0])
        t["cdf"] = np.append(t["cdf"], nxt["cdf"][0])

    cases = [("少切一张表", drop_last, "[A-1]"),
             ("cm 存错", cm_false, "[A-4]"),
             ("0.5 MeV 换成各向同性", iso_05, "[A-3]"),
             ("一张表 pdf / cdf 反", swap_rows, "[A-2]"),
             ("切片终点多 +1", shift_one, "[A-2]")]
    for what, bad, tag_want in cases:
        a = copy.deepcopy(angle)
        bad(a)
        try:
            check_A(a, quiet=True, reader=reader)
        except AssertionError as e:
            tag = str(e)[:5]
            assert tag == tag_want, "自检：「%s」响的是 %s，应是 %s" % (what, tag, tag_want)
            print("自检 ✅ %-14s → %s 响了" % (what, tag))
            continue
        raise AssertionError("自检：「%s」没有任何一条响 —— 判据有漏洞，截图给 Claude" % what)
    print("判据 A 自检全过：5 种错都响了，而且响的是对的那一条")
