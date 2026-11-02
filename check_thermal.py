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
