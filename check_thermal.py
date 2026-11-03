# -*- coding: utf-8 -*-
# ============================================================
#  check_thermal.py  --  M3 判据 T1：load_thermal 读库、切表对不对
#  Claude 2026-11-01 写（判据归 Claude，成品直接粘贴 / import，不用默写）
#
#  用法：在你写 load_thermal 的那个文件里，if __name__ == "__main__": 块中
#      from check_thermal import check_T1, check_T1_selftest
#      thermal = load_thermal(<路径>, "294K")
#      check_T1(thermal)              # 正式判
#      check_T1_selftest(thermal)     # 头一次跑做一遍：把你的结果故意改错 6 种，确认 5 条每条都会响（纪律 4）
#  运行：Windows 的 openmc (3.10.20) 环境，Ctrl+F5（只用 h5py + numpy）
#
#  两条路：
#    你的路 = load_thermal 按 energy_out 第 4 行（μ 起点）切好的 118 份小表
#    判据的路 = 这里自己直接读库，不切表；μ 起点用「全局点号 × 20」（10/31 实测两者相等）
#  门槛在看数据之前定。切对了它还响 → 截图给我，一起看是不是门槛的问题；不许直接改门槛。
# ============================================================

import os
import h5py
import numpy as np

SAB_PATH = os.path.join(r"E:\BNCT_GPU\lib80x_hdf5", "c_H_in_H2O.h5")
T_TAG = "294K"
N_MU = 20            # 每个出射点 20 个等概率 μ（10/31 实测）
E_MAX_REF = 10.0     # energy_max，eV（10/31 实测）
SIG_0253_REF = 51.69 # 0.0253 eV 的非弹截面，b（10/31 用 np.interp 读出，两位小数）


def read_raw(path=SAB_PATH, T=T_TAG):
    """判据自己的一条路：直接读库，不切表。"""
    with h5py.File(path, "r") as f:
        root = f["c_H_in_H2O"]
        dist = root[T]["inelastic/distribution"]
        raw = {
            "E_max": float(root.attrs["energy_max"]),
            "E_in": dist["energy"][:],
            "offsets": np.asarray(dist["energy_out"].attrs["offsets"]),
            "n_col": int(dist["energy_out"].shape[1]),
            "mu_row0": dist["mu"][0, :],
        }
    return raw


def check_T1(thermal, raw=None, quiet=False):
    if raw is None:
        raw = read_raw()
    tabs = thermal["tables"]
    n_tab = len(raw["E_in"])

    # [T1-1] 抓：range 少一张 / 多一张、入射能量读错数据集
    assert len(tabs) == n_tab and np.array_equal(np.asarray(thermal["E_in"]), raw["E_in"]), \
        "[T1-1] 表的张数 / 入射能量和库不一致：你 %d 张，库里 %d 张" % (len(tabs), n_tab)

    # [T1-5] 抓：xs 第 0 行 / 第 1 行拿反、energy_max 从错的地方取
    sig = np.interp(0.0253, np.asarray(thermal["xs_E"]), np.asarray(thermal["xs_sig"]))
    assert abs(thermal["E_max"] - E_MAX_REF) < 1e-9 and abs(sig - SIG_0253_REF) < 0.006, \
        "[T1-5] 截面或 energy_max 读错：E_max = %g（应 %g），σ(0.0253 eV) = %.4f b（应 %.2f）" \
        % (thermal["E_max"], E_MAX_REF, sig, SIG_0253_REF)

    worst_end = 0.0      # 最差 |cdf 末尾 − 1|
    worst_int = 0.0      # 最差 |∫pdf dE′ − (cdf 末 − cdf 首)|
    for l in range(n_tab):
        tb = tabs[l]
        E_out = np.asarray(tb["E_out"], dtype=float)
        pdf = np.asarray(tb["pdf"], dtype=float)
        cdf = np.asarray(tb["cdf"], dtype=float)
        mu = np.asarray(tb["mu"], dtype=float)
        beg = int(raw["offsets"][l])
        end = int(raw["offsets"][l + 1]) if l + 1 < n_tab else raw["n_col"]
        n = end - beg

        # [T1-2] 抓：切片起点 / 终点错一格、最后一张终点用错（cdf 会从别的表的尾巴开始，或停在半路）
        assert len(E_out) == len(pdf) == len(cdf) == n and abs(cdf[-1] - 1.0) < 1e-4 and cdf[0] < 1e-2, \
            "[T1-2] 第 %d 张切片起止不对：点数 %d（库里 %d），cdf 首 %.3g、尾 %.6f（应约 0、1）" \
            % (l, len(cdf), n, cdf[0], cdf[-1])

        # [T1-3] 抓：E′ / pdf / cdf 三行拿错（两条路：对 pdf 做梯形积分 vs 库里给的 cdf 增量）
        integ = np.trapezoid(pdf, E_out)
        assert abs(integ - (cdf[-1] - cdf[0])) < 1e-3, \
            "[T1-3] 第 %d 张 pdf 积分 %.6f ≠ cdf 增量 %.6f：E′ / pdf / cdf 三行可能拿错" \
            % (l, integ, cdf[-1] - cdf[0])

        # [T1-4] 抓：μ 起点错、μ 读成 pdf 行（全 0.05）、没拼成「点数 × 20」
        ref = raw["mu_row0"][N_MU * beg: N_MU * end].reshape(n, N_MU)
        assert mu.shape == (n, N_MU) and np.array_equal(mu, ref), \
            "[T1-4] 第 %d 张 μ 表不对：形状 %s（应 (%d, %d)），或和「全局点号 × 20」读出来的不逐位相同" \
            % (l, mu.shape, n, N_MU)

        worst_end = max(worst_end, abs(cdf[-1] - 1.0))
        worst_int = max(worst_int, abs(integ - (cdf[-1] - cdf[0])))

    if not quiet:
        print("T1 过：%d 张 · 最差 |cdf 尾 − 1| = %.2e · 最差 |∫pdf − Δcdf| = %.2e · μ 逐位相同 · "
              "σ(0.0253 eV) = %.4f b · E_max = %g eV" % (n_tab, worst_end, worst_int, sig, thermal["E_max"]))


def _with_table(thermal, l, **changes):
    """复制一份 thermal，只把第 l 张表的某几个键换掉（不动你的原件）。"""
    t = dict(thermal)
    t["tables"] = list(thermal["tables"])
    tb = dict(thermal["tables"][l])
    tb.update(changes)
    t["tables"][l] = tb
    return t


def check_T1_selftest(thermal):
    """纪律 4：在故意写错的版本上跑 T1，每种错都要响，而且响的是该响的那条。"""
    raw = read_raw()
    check_T1(thermal, raw, quiet=True)          # 原件必须先过
    tabs = thermal["tables"]
    a, b = tabs[60], tabs[61]

    def cat(x, y):
        return np.concatenate([np.asarray(x)[1:], np.asarray(y)[:1]])

    wrong = []
    t = dict(thermal); t["tables"] = list(tabs[:-1])
    wrong.append(("少了最后一张", "[T1-1]", t))
    wrong.append(("第 60 张起点右移一格", "[T1-2]",
                  _with_table(thermal, 60, E_out=cat(a["E_out"], b["E_out"]), pdf=cat(a["pdf"], b["pdf"]),
                              cdf=cat(a["cdf"], b["cdf"]),
                              mu=np.vstack([np.asarray(a["mu"])[1:], np.asarray(b["mu"])[:1]]))))
    last = tabs[-1]
    wrong.append(("最后一张终点少一列", "[T1-2]",
                  _with_table(thermal, len(tabs) - 1, E_out=np.asarray(last["E_out"])[:-1],
                              pdf=np.asarray(last["pdf"])[:-1], cdf=np.asarray(last["cdf"])[:-1],
                              mu=np.asarray(last["mu"])[:-1])))
    wrong.append(("第 30 张 pdf 用了 E′ 那一行", "[T1-3]",
                  _with_table(thermal, 30, pdf=np.asarray(tabs[30]["E_out"]).copy())))
    wrong.append(("第 10 张 μ 起点错一个出射点", "[T1-4]",
                  _with_table(thermal, 10, mu=np.vstack([np.asarray(tabs[10]["mu"])[1:],
                                                         np.asarray(tabs[11]["mu"])[:1]]))))
    t = dict(thermal); t["xs_E"], t["xs_sig"] = thermal["xs_sig"], thermal["xs_E"]
    wrong.append(("截面两行对调", "[T1-5]", t))

    for what, tag, t in wrong:
        try:
            check_T1(t, raw, quiet=True)
        except AssertionError as e:
            assert str(e).startswith(tag), "自检：「%s」响的是 %s，应该是 %s" % (what, str(e)[:6], tag)
            print("自检 ✅ %-22s → %s 响了" % (what, tag))
            continue
        raise AssertionError("自检：「%s」没有任何一条响 —— 判据有漏洞，截图给 Claude" % what)
    print("T1 自检全过：6 种错，5 条判据各自都响过")


# ============================================================
#  T3（能量那一半，11/2 Claude 加）：sample_thermal 抽出的 E′ 对不对
#  用法（thermal.py 的 if __name__ 块里）：
#      from check_thermal import check_T3E, check_T3E_selftest
#      check_T3E(thermal, sample_thermal)              # 正式判
#      check_T3E_selftest(thermal, sample_thermal)     # 头一次跑做一遍：故意写错 2 种，确认会响（纪律 4）
#  两条路：
#    你的路   = sample_thermal 抽 10 万次，数上散射比例、算 ⟨E′⟩/E
#    判据的路 = 这里直接读库里这张表的 E′、pdf 两行（不经过 load_thermal），逐段精确积分：
#               pdf 每段是直线 → ∫p 用梯形精确；E′·p 是二次 → 用辛普森精确（11/1 你推过同一个式子）
#  只在格点上测：E 取某张表自己的入射能量 → 选表 f = 0、修正 E′ 不起作用。
#      盲区：「修正 E′」那一步、两张表之间的 E —— 由整体对拍（判据 3 丙）管
#  门槛（看数据之前定）：先查 R = SE / 均值 < 0.2（规则 14），再 |z| < 3
#  明天（2）加 μ 以后，⟨μ⟩ 那一条再加进来（成为完整的 T3）
#  sample_thermal 返回一个数（E′）或一个元组（第 0 个是 E′）都能判
# ============================================================
import math
from rng import lcg

T3_N = 100000
T3_SEED = 20261103
T3_E_NEAR = [0.0253, 1.05]      # 测这两个能量附近的那张表（10/31 直接积分：54.5% / 1.19、4.7% / 0.566）
T3_Z = 3.0


def t3_read(l, path=SAB_PATH, T=T_TAG):
    """判据自己的路：直接从库里切第 l 张表的 E′、pdf 两行。"""
    with h5py.File(path, "r") as f:
        dist = f["c_H_in_H2O"][T]["inelastic/distribution"]
        E_in = dist["energy"][:]
        offs = np.asarray(dist["energy_out"].attrs["offsets"]).astype(int)
        n_col = int(dist["energy_out"].shape[1])
        beg = offs[l]
        end = offs[l + 1] if l + 1 < len(offs) else n_col
        x = dist["energy_out"][0, beg:end]
        p = dist["energy_out"][1, beg:end]
    return float(E_in[l]), np.asarray(x, float), np.asarray(p, float)


def t3_ref(E, x, p):
    """对表直接积分：上散射比例 P(E′ > E)、⟨E′⟩ / E。"""
    a, b = x[:-1], x[1:]
    pa, pb = p[:-1], p[1:]
    seg = (b - a) * (pa + pb) / 2.0                           # 每段 ∫p（梯形，对直线精确）
    norm = seg.sum()
    mean = np.sum((b - a) / 6.0 * ((2 * a + b) * pa + (a + 2 * b) * pb)) / norm
    k = int(np.searchsorted(x, E, side="right")) - 1          # E 落在第 k 段
    below = seg[:k].sum() if k > 0 else 0.0
    if 0 <= k < len(x) - 1:
        dd = E - x[k]
        s = (p[k + 1] - p[k]) / (x[k + 1] - x[k])
        below += p[k] * dd + 0.5 * s * dd * dd
    return 1.0 - below / norm, mean / E


def t3_sample(thermal, sample, E, n, seed):
    rng = lcg(seed)
    out = np.empty(n)
    for k in range(n):
        r = sample(thermal, E, rng)
        out[k] = r[0] if isinstance(r, tuple) else r
    up = float(np.mean(out > E))
    se_up = math.sqrt(up * (1 - up) / n)
    ratio = out / E
    m = float(ratio.mean())
    se_m = float(ratio.std(ddof=1)) / math.sqrt(n)
    return up, se_up, m, se_m


def check_T3E(thermal, sample, n=T3_N, quiet=False, reader=t3_read):
    E_in = thermal["E_in"]
    for t, E_near in enumerate(T3_E_NEAR):
        l = int(np.argmin(np.abs(np.asarray(E_in) - E_near)))
        E, x, p = reader(l)
        up_ref, m_ref = t3_ref(E, x, p)
        up, se_up, m, se_m = t3_sample(thermal, sample, E, n, T3_SEED + t)
        # 规则 14：先查 R
        assert se_up / up < 0.2 and se_m / m < 0.2, \
            "[T3-0] 第 %d 张 R 太大（上散射 R = %.3f，⟨E′⟩/E R = %.3f）：SE 不可信，别看 z" % (l, se_up / up, se_m / m)
        z_up = (up - up_ref) / se_up
        z_m = (m - m_ref) / se_m
        # [T3-1] 抓：选错表、区间 j 错一格、二次方程解错（E′ 偏到 E 的另一侧）
        assert abs(z_up) < T3_Z, \
            "[T3-1] 第 %d 张（E = %.4g eV）上散射比例 %.4f ± %.4f，对表积分 %.4f，z = %+.1f" % (l, E, up, se_up, up_ref, z_up)
        # [T3-2] 抓：E′ 整体偏高 / 偏低（d 少一项、修正多做了、单位）
        assert abs(z_m) < T3_Z, \
            "[T3-2] 第 %d 张（E = %.4g eV）⟨E′⟩/E = %.4f ± %.4f，对表积分 %.4f，z = %+.1f" % (l, E, m, se_m, m_ref, z_m)
        if not quiet:
            print("T3E 第 %3d 张 E = %.4g eV · 上散射 %.4f ± %.4f（表 %.4f，z = %+.2f）· ⟨E′⟩/E %.4f ± %.4f（表 %.4f，z = %+.2f）"
                  % (l, E, up, se_up, up_ref, z_up, m, se_m, m_ref, z_m))
    if not quiet:
        print("T3E 过：2 个能量 × %d 次" % n)


def check_T3E_selftest(thermal, sample, n=T3_N, reader=t3_read):
    """纪律 4：在故意写错的版本上跑，必须响。"""
    E_in = np.asarray(thermal["E_in"])

    def left_table(th, E, rng):          # 选表永远取左边一张（等价于 E 给错了一格）
        i = int(np.searchsorted(E_in, E, side="right")) - 1
        return sample(th, E_in[max(i - 1, 0)], rng)

    def high_2pct(th, E, rng):           # E′ 整体偏高 2%
        r = sample(th, E, rng)
        return (r[0] if isinstance(r, tuple) else r) * 1.02

    for what, fake in [("选表取了左边一张", left_table), ("E′ 整体偏高 2%", high_2pct)]:
        try:
            check_T3E(thermal, fake, n=n, quiet=True, reader=reader)
        except AssertionError as e:
            tag = str(e)[:6]
            assert tag in ("[T3-1]", "[T3-2]"), "自检：「%s」响的是 %s" % (what, tag)
            print("自检 ✅ %-12s → %s 响了" % (what, tag))
            continue
        raise AssertionError("自检：「%s」没有任何一条响 —— 判据有漏洞，截图给 Claude" % what)
    print("T3E 自检全过：2 种错都响了")


# ============================================================
#  T3 完整版（11/3 Claude 加）：sample_thermal 返回 (E′, μ) 以后的完整判据
#  用法（thermal.py 的 if __name__ 块里；第 13 行 import 那行末尾加上 ,check_T3,check_T3_selftest）：
#      check_T3(thermal, sample_thermal)              # 正式判（新种子 = 11/2 T3E 的「换种子复跑」）
#      check_T3_selftest(thermal, sample_thermal)     # 头一次跑做一遍：故意写错 4 种，确认每种都响（纪律 4）
#  旧的 check_T3E / check_T3E_selftest 留着不删（它们也认元组），以后用 check_T3 就行
#
#  4 条（另有 [T3-0] 前提：返回值是两个数、规则 14 的 R 预查）：
#    [T3-1] 上散射比例    [T3-2] ⟨E′⟩/E    —— 和 T3E 一样，判据的路 = 直接读库对表精确积分
#    [T3-3] ⟨μ⟩          判据的路 = 直接读库里的 μ（按「全局点号 × 20」取，不用 load_thermal 的切表、不用第 4 行）
#                         ⟨μ⟩ = Σ_j P_j · (μ̄_j + μ̄_{j+1}) / 2
#                         来历：ξ₁ 落进第 j 段的概率 P_j = 这段 pdf 的梯形积分 / 总积分；段内 ξ₁ 均匀 → g 均匀 → 平均 1/2；
#                               k 等概率 → 第 k 个 μ 平均下来 = 20 个的平均 μ̄；抹开左右对称 → 不改平均
#    [T3-4] 算例逐位     假发生器喂固定的 ξ₁、ξ₂、ξ₃：第 48 张、第 300 个出射点，
#                         k = 3（中间）· 0（左边镜像）· 19（右边镜像）各一次（g = 0），再 k = 3、g ≈ 0.37 一次；
#                         和这里直接从库算的逐位比（容差 1e-12）；顺带查一次热散射正好用 3 个随机数
#  门槛（看数据之前定）：T3-1～3 |z| < 3；T3-4 容差 1e-12
#  盲区：两张表之间的 E、修正 E′ —— 整体对拍（判据 3 丙）管；g 写成 1 − g 时 T3-3 抓不到（平均还是 1/2），靠 T3-4 的 g ≈ 0.37 那次
#  ⟨μ⟩ 不做 R 预查：μ 有正有负，均值接近 0 时 R = SE / 均值 没有意义；10 万次的 SE 本身可信
# ============================================================

T3_SEED_NEW = 20261104          # 11/3 换种子（11/2 的 T3E 用 20261103）
T3_CASE_L = 48                  # 讲义 §5.5 算例：第 48 张表（0.0253 eV）
T3_CASE_J = 300                 #                第 300 个出射点
T3_CASE = [                     # (g 的目标值, ξ₂, ξ₃)
    (0.0, 0.17, 0.73),          # k = 3：讲义 §5.5 那个算例，讲义写 −0.194259
    (0.0, 0.01, 0.73),          # k = 0：左边没有邻居 → 镜像
    (0.0, 0.98, 0.73),          # k = 19：右边没有邻居 → 镜像
    (0.37, 0.17, 0.73),         # k = 3，g ≈ 0.37：查 g 的方向（g vs 1 − g）
]
T3_TOL = 1e-12


class FakeRng:
    """假发生器：按顺序吐出给定的数；多要一个就响（一次热散射只该用 ξ₁、ξ₂、ξ₃ 三个）。"""
    def __init__(self, seq):
        self.seq = list(seq)
        self.used = 0

    def random(self):
        assert self.used < len(self.seq), \
            "[T3-4] sample_thermal 要了第 %d 个随机数：一次热散射只该用 3 个（ξ₁ 是不是又抽了一次？）" % (self.used + 1)
        x = self.seq[self.used]
        self.used += 1
        return x


def t3_read_full(l, path=SAB_PATH, T=T_TAG):
    """判据自己的路：第 l 张表的 E′、pdf、cdf 三行 + μ（按全局点号 × 20 取，形状 点数 × 20）。"""
    with h5py.File(path, "r") as f:
        dist = f["c_H_in_H2O"][T]["inelastic/distribution"]
        E_in = dist["energy"][:]
        offs = np.asarray(dist["energy_out"].attrs["offsets"]).astype(int)
        n_col = int(dist["energy_out"].shape[1])
        beg = int(offs[l])
        end = int(offs[l + 1]) if l + 1 < len(offs) else n_col
        x = dist["energy_out"][0, beg:end]
        p = dist["energy_out"][1, beg:end]
        c = dist["energy_out"][2, beg:end]
        mu = dist["mu"][0, beg * N_MU:end * N_MU]
    return (float(E_in[l]), np.asarray(x, float), np.asarray(p, float),
            np.asarray(c, float), np.asarray(mu, float).reshape(-1, N_MU))


def t3_mu_ref(x, p, mu):
    """对表直接算 ⟨μ⟩ = Σ_j P_j · (μ̄_j + μ̄_{j+1}) / 2。"""
    seg = (x[1:] - x[:-1]) * (p[:-1] + p[1:]) / 2.0
    mbar = mu.mean(axis=1)
    return float(np.sum(seg * (mbar[:-1] + mbar[1:]) / 2.0) / seg.sum())


def t3_case_want(c, mu, j, xi1, xi2, xi3):
    """算例的「应该得多少」：照讲义 §5.5 一步一步算（判据自己写的一份）。"""
    g = (xi1 - c[j]) / (c[j + 1] - c[j])
    k = int(xi2 * N_MU)

    def at(kk):
        return mu[j, kk] + g * (mu[j + 1, kk] - mu[j, kk])

    m = at(k)
    left = at(k - 1) if k > 0 else -1.0 - (m + 1.0)
    right = at(k + 1) if k < N_MU - 1 else 1.0 + (1.0 - m)
    return m + min(m - left, right - m) * (xi3 - 0.5), k, g


def check_T3(thermal, sample, n=T3_N, seed=T3_SEED_NEW, quiet=False, reader=t3_read_full):
    E_in = np.asarray(thermal["E_in"])

    # ---------- 统计部分：两张表各抽 n 次 ----------
    for t, E_near in enumerate(T3_E_NEAR):
        l = int(np.argmin(np.abs(E_in - E_near)))
        E, x, p, c, mu = reader(l)
        up_ref, m_ref = t3_ref(E, x, p)
        mu_ref = t3_mu_ref(x, p, mu)
        rng = lcg(seed + t)
        Eo = np.empty(n)
        mo = np.empty(n)
        for i in range(n):
            r = sample(thermal, E, rng)
            assert isinstance(r, tuple) and len(r) == 2, \
                "[T3-0] sample_thermal 应返回 (E′, μ) 两个数，现在返回的是 %r" % (r,)
            Eo[i], mo[i] = r
        up = float(np.mean(Eo > E))
        se_up = math.sqrt(up * (1 - up) / n)
        ratio = Eo / E
        m = float(ratio.mean())
        se_m = float(ratio.std(ddof=1)) / math.sqrt(n)
        mu_mean = float(mo.mean())
        se_mu = float(mo.std(ddof=1)) / math.sqrt(n)
        # 规则 14：先查 R（只对上散射比例、⟨E′⟩/E）
        assert se_up / up < 0.2 and se_m / m < 0.2, \
            "[T3-0] 第 %d 张 R 太大（上散射 R = %.3f，⟨E′⟩/E R = %.3f）：SE 不可信，别看 z" % (l, se_up / up, se_m / m)
        z_up = (up - up_ref) / se_up
        z_m = (m - m_ref) / se_m
        z_mu = (mu_mean - mu_ref) / se_mu
        # [T3-1] 抓：选错表、区间 j 错一格、二次方程解错（E′ 偏到 E 的另一侧）
        assert abs(z_up) < T3_Z, \
            "[T3-1] 第 %d 张（E = %.4g eV）上散射比例 %.4f ± %.4f，对表积分 %.4f，z = %+.1f" % (l, E, up, se_up, up_ref, z_up)
        # [T3-2] 抓：E′ 整体偏高 / 偏低（d 少一项、修正多做了）
        assert abs(z_m) < T3_Z, \
            "[T3-2] 第 %d 张（E = %.4g eV）⟨E′⟩/E = %.4f ± %.4f，对表积分 %.4f，z = %+.1f" % (l, E, m, se_m, m_ref, z_m)
        # [T3-3] 抓：μ 符号反、取错行 / 列、μ 表和出射点对不上、抹开不对称
        assert abs(z_mu) < T3_Z, \
            "[T3-3] 第 %d 张（E = %.4g eV）⟨μ⟩ = %.4f ± %.4f，对表直接算 %.4f，z = %+.1f" % (l, E, mu_mean, se_mu, mu_ref, z_mu)
        if not quiet:
            print("T3 第 %3d 张 E = %.4g eV · 上散射 %.4f ± %.4f（表 %.4f，z = %+.2f）· ⟨E′⟩/E %.4f ± %.4f（表 %.4f，z = %+.2f）"
                  % (l, E, up, se_up, up_ref, z_up, m, se_m, m_ref, z_m))
            print("               · ⟨μ⟩ %.4f ± %.4f（表 %.4f，z = %+.2f）· μ 范围 [%.5f, %.5f]"
                  % (mu_mean, se_mu, mu_ref, z_mu, mo.min(), mo.max()))

    # ---------- 算例逐位：假发生器 ----------
    E, x, p, c, mu = reader(T3_CASE_L)
    j = T3_CASE_J
    assert c[j + 1] > c[j], "判据自己的问题：第 %d 张第 %d 段 cdf 不增，换一个点（截图给 Claude）" % (T3_CASE_L, j)
    for gt, xi2, xi3 in T3_CASE:
        xi1 = c[j] + gt * (c[j + 1] - c[j])
        want, k, g = t3_case_want(c, mu, j, xi1, xi2, xi3)
        fake = FakeRng([xi1, xi2, xi3])
        r = sample(thermal, E, fake)
        got = r[1]
        # [T3-4] 抓：k 取整错、[j, k] 行列反、k = 0 / 19 没镜像（k − 1 = −1 悄悄取到第 19 个）、
        #        g 用错（1 − g、按能量插）、ξ₂ / ξ₃ 用反、ξ₁ 又抽了一次
        assert fake.used == 3, \
            "[T3-4] 一次热散射只用了 %d 个随机数（应 3 个：ξ₁、ξ₂、ξ₃）" % fake.used
        assert abs(got - want) < T3_TOL, \
            "[T3-4] 算例 k = %d、g = %.2f：你的 μ = %.10f，应为 %.10f（差 %.2e）" % (k, g, got, want, got - want)
        if not quiet:
            print("T3 算例 第 %d 张第 %d 点 · k = %2d · g = %.2f · μ = %+.7f（应 %+.7f）✅" % (T3_CASE_L, j, k, g, got, want))
    if not quiet:
        print("T3 过：2 个能量 × %d 次（种子 %d、%d）+ 4 个算例逐位（讲义 §5.5 写 −0.194259，对第一个算例）"
              % (n, seed, seed + 1))


def check_T3_selftest(thermal, sample, n=T3_N, reader=t3_read_full):
    """纪律 4：在故意写错的版本上跑，必须响，而且要响对的那一条。"""
    E_in = np.asarray(thermal["E_in"])

    def left_table(th, E, rng):          # 选表永远取左边一张 → 应 T3-1 / T3-2 响
        i = int(np.searchsorted(E_in, E, side="right")) - 1
        return sample(th, E_in[max(i - 1, 0)], rng)

    def high_2pct(th, E, rng):           # E′ 整体偏高 2% → 应 T3-1 / T3-2 响
        Ep, mu = sample(th, E, rng)
        return Ep * 1.02, mu

    def mu_flip(th, E, rng):             # μ 符号反了 → 应 T3-3 响
        Ep, mu = sample(th, E, rng)
        return Ep, -mu

    def mu_nudge(th, E, rng):            # μ 只差 1e-7（统计上看不出）→ 应 T3-4 响
        Ep, mu = sample(th, E, rng)
        return Ep, mu + 1e-7

    cases = [("选表取了左边一张", left_table, ("[T3-1]", "[T3-2]")),
             ("E′ 整体偏高 2%", high_2pct, ("[T3-1]", "[T3-2]")),
             ("μ 符号反了", mu_flip, ("[T3-3]",)),
             ("μ 差 1e-7", mu_nudge, ("[T3-4]",))]
    for what, fake, tags in cases:
        try:
            check_T3(thermal, fake, n=n, quiet=True, reader=reader)
        except AssertionError as e:
            tag = str(e)[:6]
            assert tag in tags, "自检：「%s」响的是 %s，应是 %s" % (what, tag, " / ".join(tags))
            print("自检 ✅ %-12s → %s 响了" % (what, tag))
            continue
        raise AssertionError("自检：「%s」没有任何一条响 —— 判据有漏洞，截图给 Claude" % what)
    print("T3 自检全过：4 种错都响了，而且响的是对的那一条")


# ============================================================
#  T2（11/3 Claude 加）：rotate_dir(u, v, w, mu, phi) 转方向对不对
#  用法（thermal.py 的 if __name__ 块里；import 那行再加上 ,check_T2,check_T2_selftest）：
#      check_T2(rotate_dir)              # 正式判
#      check_T2_selftest()               # 头一次跑做一遍：判据在 5 个故意写错的版本上必须响（纪律 4）
#  不读库、不抽随机数，几秒钟跑完
#
#  两条路：
#    你的路   = rotate_dir 逐分量的三行公式（讲义 §5.6 / OpenMC rotate_angle）
#    判据的路 = 向量写法：新方向 = μ · 旧方向 + a · (cos φ · e₁ + sin φ · e₂)，
#               e₁、e₂ 是和旧方向垂直的两根单位向量（讲义 §5.6「来历」那两根；特殊分支换成绕 v 的那两根），
#               用 numpy 整体算，不逐分量抄公式
#  4 条（跑的顺序 T2-2 → T2-3 → T2-1 → T2-4）：
#    [T2-1] 逐位：讲义 §5.6 算例 + 6 组一般方向，和判据的路比（容差 1e-12）
#    [T2-2] 长度 = 1：2000 组随机 (u, v, w, μ, φ)，|新方向|² − 1 < 1e-9
#    [T2-3] 夹角 = μ：同 2000 组，旧 · 新 − μ < 1e-9
#    [T2-4] w = ±1 不出 nan：旧方向 (0, 0, ±1) 和 (很小, 很小, ±1)，结果有限、和判据的路逐位（走特殊分支）
#  门槛：逐位 1e-12；长度、夹角 1e-9（11/3 Claude 用照 OpenMC 写的对照版试跑：随机方向里 w = 0.99999 那种，
#        1 − w² 掉有效位，长度误差到 2e-11 —— 门槛是在看到这个之后从 1e-12 放到 1e-9 的，照实记）
#  盲区：φ 的方向约定（φ → −φ 时 T2-2、T2-3 抓不到）只有 T2-1 管 —— 所以 T2-1 用了 7 组
# ============================================================

T2_TOL = 1e-12       # 逐位（T2-1、T2-4）
T2_TOL_GEO = 1e-9    # 随机方向的长度、夹角（T2-2、T2-3）：w 接近 ±1 时 1 − w² 掉有效位，实测最差约 2e-11
T2_N = 2000
T2_SEED = 20261103


def t2_want(u, v, w, mu, phi):
    """判据的路：向量写法（不逐分量抄公式）。"""
    d = np.array([u, v, w], float)
    a = math.sqrt(max(0.0, 1.0 - mu * mu))
    b = math.sqrt(max(0.0, 1.0 - w * w))
    if b > 1e-10:
        e1 = np.array([u * w, v * w, -(1.0 - w * w)]) / b     # 讲义 §5.6 来历里的 e₁
        e2 = np.array([-v, u, 0.0]) / b                       #                    e₂
    else:                                                     # 绕 v 展开（OpenMC rotate_angle 的特殊分支）
        b = math.sqrt(1.0 - v * v)
        e1 = np.array([u * v, -(1.0 - v * v), v * w]) / b
        e2 = np.array([w, 0.0, -u]) / b
    return mu * d + a * (math.cos(phi) * e1 + math.sin(phi) * e2)


def t2_cases_exact():
    """T2-1 的 7 组：第 0 组 = 讲义 §5.6 算例；其余 6 组各方向都斜、μ、φ 都不取特殊值。"""
    cases = [(0.36, 0.48, 0.80, -0.194259, 2.1)]
    raw = [(0.27, -0.54, 0.31, 0.43, 0.77), (-0.62, 0.18, -0.49, -0.71, 4.4),
           (0.11, 0.93, -0.22, 0.86, 3.3), (-0.38, -0.47, 0.66, -0.29, 5.9),
           (0.71, 0.08, -0.53, 0.12, 1.4), (-0.24, 0.39, 0.83, 0.57, 2.8)]
    for x, y, z, mu, phi in raw:
        n = math.sqrt(x * x + y * y + z * z)
        cases.append((x / n, y / n, z / n, mu, phi))
    return cases


def t2_cases_random(n=T2_N, seed=T2_SEED):
    rs = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        d = rs.normal(size=3)
        d /= np.linalg.norm(d)
        out.append((float(d[0]), float(d[1]), float(d[2]), float(rs.uniform(-1, 1)), float(rs.uniform(0, 2 * math.pi))))
    return out


def t2_cases_pole():
    """旧方向贴着 ±z 轴：正好 (0, 0, ±1)，以及 u、v 不为 0 但 w 仍是 ±1.0（1 − w² 算出来正好 0 → 特殊分支）。
    注：双精度里 w 离 1 最近也差 1.1e-16 → 1 − w² ≥ 2.2e-16 → b ≥ 1.5e-8，所以 b 落在 (0, 1e-10] 的情况造不出来。"""
    out = []
    for s in (1.0, -1.0):
        out.append((0.0, 0.0, s, 0.38, 1.9))
        out.append((0.0, 0.0, s, -0.64, 4.7))
        out.append((2e-9, -3e-9, s, -0.17, 3.6))
    return out


def check_T2(rotate, quiet=False):
    worst_len = 0.0
    worst_dot = 0.0
    for u, v, w, mu, phi in t2_cases_random():
        un, vn, wn = rotate(u, v, w, mu, phi)
        dl = abs(un * un + vn * vn + wn * wn - 1.0)
        dd = abs(u * un + v * vn + w * wn - mu)
        worst_len = max(worst_len, dl)
        worst_dot = max(worst_dot, dd)
        # [T2-2] 抓：少除 b、a 写成 √(1 − w²)、某一项漏乘 a
        assert dl < T2_TOL_GEO, "[T2-2] 旧方向 (%.4f, %.4f, %.4f)、μ = %.4f、φ = %.3f：新方向长度² = %.15f ≠ 1" % (u, v, w, mu, phi, dl + 1.0)
        # [T2-3] 抓：μ 那一项漏写 / 乘错分量、新方向和旧方向的夹角不是 μ
        assert dd < T2_TOL_GEO, "[T2-3] 旧方向 (%.4f, %.4f, %.4f)、μ = %.4f、φ = %.3f：旧 · 新 = %.12f ≠ μ" % (u, v, w, mu, phi, dd + mu)
    if not quiet:
        print("T2-2、T2-3 ✅ %d 组随机方向 · 最差 |长度² − 1| = %.1e · 最差 |旧·新 − μ| = %.1e" % (T2_N, worst_len, worst_dot))

    # [T2-1] 抓：抄错一个分量、正负号错（包括 φ → −φ 这种长度、夹角都对的错）、u / v 写反
    for u, v, w, mu, phi in t2_cases_exact():
        got = np.array(rotate(u, v, w, mu, phi), float)
        want = t2_want(u, v, w, mu, phi)
        assert np.all(np.abs(got - want) < T2_TOL), \
            "[T2-1] 旧方向 (%.4f, %.4f, %.4f)、μ = %.4f、φ = %.2f：你的 (%.6f, %.6f, %.6f)，应 (%.6f, %.6f, %.6f)" \
            % (u, v, w, mu, phi, got[0], got[1], got[2], want[0], want[1], want[2])
    if not quiet:
        u, v, w, mu, phi = t2_cases_exact()[0]
        g = rotate(u, v, w, mu, phi)
        print("T2-1 ✅ 7 组逐位 · 讲义 §5.6 算例：(%.5f, %.5f, %.5f)（讲义 (−0.98506, 0.09787, 0.14173)）" % tuple(g))

    for u, v, w, mu, phi in t2_cases_pole():
        got = np.array(rotate(u, v, w, mu, phi), float)
        want = t2_want(u, v, w, mu, phi)
        # [T2-4] 抓：没写特殊分支（除以 0 → nan / ZeroDivisionError）、特殊分支公式抄错
        assert np.all(np.isfinite(got)), "[T2-4] 旧方向 (%.1e, %.1e, %+.1f)：结果有 nan / inf %r —— 特殊分支没进去？" % (u, v, w, tuple(got))
        assert np.all(np.abs(got - want) < T2_TOL), \
            "[T2-4] 旧方向 (%.1e, %.1e, %+.1f)、μ = %.2f、φ = %.1f：你的 (%.6f, %.6f, %.6f)，应 (%.6f, %.6f, %.6f)" \
            % (u, v, w, mu, phi, got[0], got[1], got[2], want[0], want[1], want[2])
    if not quiet:
        print("T2-4 ✅ 旧方向贴着 ±z 轴的 6 组：没有 nan、和判据的路一致")
        print("T2 过")


def check_T2_selftest():
    """纪律 4：判据在故意写错的 rotate 上必须响，而且响对的那一条。"""
    def good(u, v, w, mu, phi):                   # 判据自己的对照：向量写法
        return tuple(t2_want(u, v, w, mu, phi))

    def phi_flip(u, v, w, mu, phi):               # φ → −φ：长度、夹角都对 → 只有 T2-1 能抓
        return good(u, v, w, mu, -phi)

    def no_div_b(u, v, w, mu, phi):               # 一般分支少除 b
        a = math.sqrt(1 - mu * mu)
        b = math.sqrt(max(0.0, 1 - w * w))
        c, s = math.cos(phi), math.sin(phi)
        return (mu * u + a * (u * w * c - v * s), mu * v + a * (v * w * c + u * s), mu * w - a * b * c)

    def mu_sign(u, v, w, mu, phi):                # μ 那一项符号反了：长度还是 1，夹角变成 −μ
        g = t2_want(u, v, w, mu, phi)
        return tuple(g - 2 * mu * np.array([u, v, w]))

    def no_pole(u, v, w, mu, phi):                # 没写特殊分支
        a = math.sqrt(1 - mu * mu)
        b = math.sqrt(max(0.0, 1 - w * w))
        c, s = math.cos(phi), math.sin(phi)
        with np.errstate(all="ignore"):
            b = np.float64(b)
            return (mu * u + a * (u * w * c - v * s) / b, mu * v + a * (v * w * c + u * s) / b, mu * w - a * b * c)

    def pole_typo(u, v, w, mu, phi):              # 特殊分支 w′ 那一行 − u sin φ 写成 + u sin φ
        b = math.sqrt(max(0.0, 1 - w * w))
        if b > 1e-10:
            return good(u, v, w, mu, phi)
        a = math.sqrt(1 - mu * mu)
        b = math.sqrt(1 - v * v)
        c, s = math.cos(phi), math.sin(phi)
        return (mu * u + a * (u * v * c + w * s) / b, mu * v - a * b * c, mu * w + a * (v * w * c + u * s) / b)

    check_T2(good, quiet=True)                    # 对照本身必须过，否则判据自己有问题
    cases = [("φ 方向反了", phi_flip, ("[T2-1]",)),
             ("少除 b", no_div_b, ("[T2-2]",)),
             ("μ 那一项符号反", mu_sign, ("[T2-3]",)),
             ("没写特殊分支", no_pole, ("[T2-4]",)),
             ("特殊分支抄错号", pole_typo, ("[T2-4]",))]
    for what, fake, tags in cases:
        try:
            check_T2(fake, quiet=True)
        except AssertionError as e:
            tag = str(e)[:6]
            assert tag in tags, "自检：「%s」响的是 %s，应是 %s" % (what, tag, " / ".join(tags))
            print("自检 ✅ %-10s → %s 响了" % (what, tag))
            continue
        raise AssertionError("自检：「%s」没有任何一条响 —— 判据有漏洞，截图给 Claude" % what)
    print("T2 自检全过：5 种错都响了，而且响的是对的那一条")
