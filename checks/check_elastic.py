# -*- coding: utf-8 -*-
# ============================================================
#  check_elastic.py  --  M3 判据 A：load_elastic_angle 读库、切表对不对 · 判据 B（文件末尾）：sample_elastic_mu 抽样
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

# ---- 路径（2026-11-10 整理文件夹，Claude 加）：本文件在 my_code 的子文件夹里。
#      下面几行把 my_code 根目录和各子文件夹加进搜索路径，后面的 import 才找得到别的文件。----
import os, sys
_MY_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # my_code 根目录
for _sub in ("whiteboard", "lessons", "openmc_compare", "checks", ""):
    _p = os.path.join(_MY_CODE, _sub) if _sub else _MY_CODE
    if _p not in sys.path:
        sys.path.insert(0, _p)

import os
import copy
import math
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


# ============================================================
#  判据 B：sample_elastic_mu 抽出来的 μ_cm 分布对不对
#  Claude 2026-11-05 写（11/4 定的 3 条；11/4 的草稿没存下来，11/5 重写，种子 / N 在跑之前重定）
#
#  用法：在 elastic_angle.py 的 if __name__ == "__main__": 块里
#      from check_elastic import check_B, check_B_selftest
#      check_B(angle, sample_elastic_mu)     # 正式判
#      check_B_selftest(angle)               # 头一次跑做一遍：判据自带的正确版要过，4 种错要响对
#
#  两条路：
#    你的路   = 抽 N 次 μ_cm 求平均（随机）
#    判据的路 = 不抽样：对表精确积分得 μ̄（a_mubar，辛普森），选表按 (1 − f) : f 两张表混合
#  B-3、B-1 比的是同一个量 ⟨μ_cm⟩：先查 R = SE / |均值| < 0.2（规则 14），再看 |z| < 3
#  B-2 比的是同一个量「μ 落进每个箱的概率」：判据正着积分 pdf，你的函数反着抽 —— 两条路
#  3 条（跑的顺序 B-2 → B-3 → B-1）：
#    [B-2] 玩具三段 U 形表 μ [−1, −0.4, 0.4, 1]、pdf [0.6, 0.1, 0.3, 1.8]、cdf [0, 0.21, 0.37, 1]：
#          抽 N 次，20 个等宽箱的直方图 vs 正着积分出来的每箱概率，χ²（19 个自由度）< 43.8（p = 0.001）
#          抓：段内反解错（当均匀 / 公式错）、在 mu 里找 j（这种错格点上的 cdf 全对、平均值也几乎对，只有段内形状错）
#    [B-3] 玩具两张表，左前倾（0.1 → 0.9）右后倾（0.9 → 0.1），E 在 30% 处：
#          ⟨μ⟩ = 0.7 × 0.26667 + 0.3 × (−0.26667) = +0.10667                          抓：选表方向反、只用最近的表
#    [B-1] 真库 O16 0.5 MeV、2 MeV：⟨μ⟩ vs 表积分（11/3 判据 A：+0.4166、+0.1211）   抓：传错字典 / 核素、真表上才出现的问题
#          （自检里没有只让 B-1 响的错：B-1 是真库上的整体核对，照实记）
#  11/5 改：B-2 从 11/4 定的「单段表比平均值」换成「三段表比直方图」—— 试跑发现「在 mu 里找 j」
#          在单段表上和对的完全一样、在三段表上平均值也只差 0.001，只有看形状才抓得到
#  门槛在看数据之前定：B-3、B-1 |z| < 3、R < 0.2；B-2 χ² < 43.8；种子 B_SEED = 20261105，N = 100000（每条）
#  盲区：B-3、B-1 只看平均值；箱宽 0.1 以内的形状错 B-2 也看不到
# ============================================================

B_SEED = 20261105
B_N = 100000
B_NBIN = 20
B_CHI2_MAX = 43.8          # χ²(19) 的 99.9% 分位
_B_MU = [-1.0, 1.0]
_B_FWD = {"mu": _B_MU, "pdf": [0.1, 0.9], "cdf": [0.0, 1.0]}     # 前倾：⟨μ⟩ = +0.26667
_B_BWD = {"mu": _B_MU, "pdf": [0.9, 0.1], "cdf": [0.0, 1.0]}     # 后倾：⟨μ⟩ = −0.26667
_B_U3 = {"mu": [-1.0, -0.4, 0.4, 1.0], "pdf": [0.6, 0.1, 0.3, 1.8], "cdf": [0.0, 0.21, 0.37, 1.0]}   # U 形三段表


def b_expect(ang, E):
    """判据的路：不抽样。按 E 找左右两张表、f，⟨μ⟩ = (1 − f)·μ̄_左 + f·μ̄_右（μ̄ 用 a_mubar 精确积分）。"""
    Eg = np.asarray(ang["E"], float)
    i = int(np.clip(np.searchsorted(Eg, E, side="right") - 1, 0, len(Eg) - 2))
    f = float(np.clip((E - Eg[i]) / (Eg[i + 1] - Eg[i]), 0.0, 1.0))
    m = [a_mubar(np.asarray(ang["tables"][k]["mu"], float), np.asarray(ang["tables"][k]["pdf"], float)) for k in (i, i + 1)]
    return (1.0 - f) * m[0] + f * m[1], i, f


def b_one(tag, sampler, ang, E, rng, quiet):
    want, i, f = b_expect(ang, E)
    x = np.array([sampler(ang, E, rng) for _ in range(B_N)])
    mean = float(x.mean())
    se = float(x.std(ddof=1) / np.sqrt(B_N))
    R = se / abs(mean) if mean != 0.0 else float("inf")
    z = (mean - want) / se
    # 先查 R：样本均值离 0 太近（R > 0.2）时 SE 不可信，z 没有意义 —— 这本身就说明分布不对
    assert R < 0.2, "%s ⟨μ⟩ = %+.5f ± %.5f，R = %.2f > 0.2（应是 %+.5f）—— 平均值几乎是 0，分布不对" % (tag, mean, se, R, want)
    # z：抽样平均和精确积分差几个标准误
    assert abs(z) < 3.0, "%s ⟨μ⟩ = %+.5f ± %.5f，表积分 %+.5f，z = %+.2f（|z| 应 < 3）" % (tag, mean, se, want, z)
    if not quiet:
        print("%s ✅ E = %.4g eV（第 %d、%d 张，f = %.3f）· ⟨μ⟩ = %+.5f ± %.5f · 表积分 %+.5f · R = %.3f · z = %+.2f"
              % (tag, E, i, i + 1, f, mean, se, want, R, z))


def b_cdf_exact(t, x):
    """判据的路（正着算）：pdf 每段是直线 → F(x) = c_k + p_k·d + s·d²/2，d = x − μ_k。"""
    mu, p, c = (np.asarray(t[k], float) for k in ("mu", "pdf", "cdf"))
    k = int(np.clip(np.searchsorted(mu, x, side="right") - 1, 0, len(mu) - 2))
    s = (p[k + 1] - p[k]) / (mu[k + 1] - mu[k])
    d = x - mu[k]
    return c[k] + p[k] * d + 0.5 * s * d * d


def b_hist(tag, sampler, ang, E, rng, quiet):
    t = ang["tables"][0]
    edges = np.linspace(-1.0, 1.0, B_NBIN + 1)
    F = np.array([b_cdf_exact(t, x) for x in edges])
    expect = B_N * np.diff(F) / (F[-1] - F[0])
    x = np.array([sampler(ang, E, rng) for _ in range(B_N)])
    obs = np.histogram(x, bins=edges)[0]
    chi2 = float(np.sum((obs - expect) ** 2 / expect))
    # χ²：20 个箱里抽到的个数和正着积分的期望差多少（每箱期望都 > 5 才能用 χ²）
    assert expect.min() > 5.0, "%s 有的箱期望个数 < 5（%.1f），χ² 不能用 —— 判据本身的问题，截图给 Claude" % (tag, expect.min())
    assert chi2 < B_CHI2_MAX, "%s 直方图和表的形状对不上：χ² = %.1f（应 < %.1f，%d 个自由度）" % (tag, chi2, B_CHI2_MAX, B_NBIN - 1)
    if not quiet:
        print("%s ✅ 玩具三段 U 形表 · 20 箱直方图 χ² = %.1f（< %.1f，19 个自由度）" % (tag, chi2, B_CHI2_MAX))


def check_B(angle, sampler, quiet=False):
    from rng import lcg
    rng = lcg(B_SEED)
    toy2 = {"E": [1.0e5, 2.0e5], "cm": True, "tables": [_B_U3, _B_U3]}
    b_hist("[B-2]", sampler, toy2, 1.3e5, rng, quiet)
    toy3 = {"E": [1.0e5, 2.0e5], "cm": True, "tables": [_B_FWD, _B_BWD]}
    b_one("[B-3]", sampler, toy3, 1.3e5, rng, quiet)
    for E in (5.0e5, 2.0e6):
        b_one("[B-1]", sampler, angle["O16"], E, rng, quiet)
    if not quiet:
        print("判据 B 过")


def _b_core(ang, E, rng, pick, find_in, uniform):
    """判据自带的抽样（自检用），三个开关各对应一种错。pick: "rand" / "rev" / "near"。"""
    Eg = ang["E"]
    i = int(np.clip(np.searchsorted(Eg, E, side="right") - 1, 0, len(Eg) - 2))
    f = min(max(0.0, (E - Eg[i]) / (Eg[i + 1] - Eg[i])), 1.0)
    xi1 = rng.random()
    if pick == "rand":
        l = i + 1 if f > xi1 else i
    elif pick == "rev":
        l = i if f > xi1 else i + 1
    else:
        l = i + 1 if f > 0.5 else i
    t = ang["tables"][l]
    mu, p, c = t["mu"], t["pdf"], t["cdf"]
    xi2 = rng.random()
    j = int(np.searchsorted(t[find_in], xi2, side="right") - 1)
    j = min(max(0, j), len(c) - 2)
    if uniform:
        x = mu[j] + (xi2 - c[j]) / (c[j + 1] - c[j]) * (mu[j + 1] - mu[j])
    else:
        s = (p[j + 1] - p[j]) / (mu[j + 1] - mu[j])
        if s == 0.0:
            d = (xi2 - c[j]) / p[j]
        else:
            d = (math.sqrt(max(0.0, p[j] ** 2 + 2 * s * (xi2 - c[j]))) - p[j]) / s
        x = mu[j] + d
    return min(max(-1.0, x), 1.0)


def check_B_selftest(angle):
    """纪律 4：判据自带的正确版必须过；4 种故意写错的必须响，而且响对的那一条。"""
    ok = lambda a, E, r: _b_core(a, E, r, "rand", "cdf", False)
    check_B(angle, ok, quiet=True)
    print("自检 ✅ 判据自带的正确版 → 全过")
    cases = [("段内当均匀", lambda a, E, r: _b_core(a, E, r, "rand", "cdf", True), "[B-2]"),
             ("选表方向反", lambda a, E, r: _b_core(a, E, r, "rev", "cdf", False), "[B-3]"),
             ("只用最近的表", lambda a, E, r: _b_core(a, E, r, "near", "cdf", False), "[B-3]"),
             ("在 mu 里找 j", lambda a, E, r: _b_core(a, E, r, "rand", "mu", False), "[B-2]")]
    for what, bad, tag_want in cases:
        try:
            check_B(angle, bad, quiet=True)
        except AssertionError as e:
            tag = str(e)[:5]
            assert tag == tag_want, "自检：「%s」响的是 %s，应是 %s" % (what, tag, tag_want)
            print("自检 ✅ %-10s → %s 响了" % (what, tag))
            continue
        raise AssertionError("自检：「%s」没有任何一条响 —— 判据有漏洞，截图给 Claude" % what)
    print("判据 B 自检全过：正确版过，4 种错都响了，而且响的是对的那一条")
