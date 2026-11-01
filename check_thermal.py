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
